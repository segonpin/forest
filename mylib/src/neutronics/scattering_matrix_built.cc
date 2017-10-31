/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_matrix_built.cc
 * @brief  Implementation of class template ScatteringMatrixBuilt
 */

#include "neutronics/scattering_matrix_built.h"

#include "neutronics/scattering_source.h"
#include "neutronics/state.h"
#include "input/input.h"
#include "input/input_mat.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/base/quadrature_lib.h"  // for QGauss
#include "deal.II/fe/fe_values.h"
#include "deal.II/lac/full_matrix.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/meshworker/dof_info.h"
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/simple.h"
#include "deal.II/meshworker/loop.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/base/std_cxx11/bind.h"
#include "deal.II/dofs/dof_tools.h"     // for FiniteElement, Mapping, etc
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc
#include "deal.II/lac/block_sparse_matrix.h"  // for BlockSparseMatrix
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/lac/block_vector_base.h"  // for BlockVectorBase
#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/vector.h"         // for Vector
#include "deal.II/meshworker/local_integrator.h"  // for DoFInfo
#include "deal.II/meshworker/vector_selector.h"  // for FEValuesBase

#include <vector>                       // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  ScatteringMatrixBuilt<dim>::ScatteringMatrixBuilt (State<dim> & state)
      : mp_state (state),
        mp_fe (state.get_fe ()),
        mp_mapping (state.get_mapping ()),
        mp_dof_handler (state.get_dof_handler ())
  {
    // We first setup the non zero pattern, in order to avoid operations later
    this->set_non_zero_pattern ();

    setup_system ();
    build_system ();
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::set_non_zero_pattern ()
  {
    // we consider zero anything below this value;
    const double eps = 1.e-14;
    const unsigned int n_groups = mp_state.get_n_groups ();
    m_non_zero.resize (n_groups, std::vector<bool> (n_groups, false));
    double sigmas;
    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        for (cell = mp_dof_handler.begin_active (); cell != endc; ++cell)
        {
          /*
           * Here the sigmas should come from TH
           */
          sigmas =
              mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_sigmas[g][h];
          if (std::abs (sigmas) > eps)
          {
            m_non_zero[g][h] = true;
            // Once we find a non-zero value, then stop the loop for this group
            continue;
          }
        }
      }
    }
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::print_non_zero_pattern ()
  {
    for (unsigned int g = 0; g < m_non_zero.size (); ++g)
    {
      for (unsigned int h = 0; h < m_non_zero[g].size (); ++h)
      {
        std::cout << m_non_zero[g][h] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  template <int dim>
  bool
  ScatteringMatrixBuilt<dim>::no_upscatt (const unsigned int g) const
  {
    for (unsigned int h = g + 1; h < m_non_zero[g].size (); ++h)
    {
      if (m_non_zero[g][h])
      {
        return false;
      }
    }
    return true;
  }

  template <int dim>
  double
  ScatteringMatrixBuilt<dim>::memory_consumption () const
  {
    double mem_consumtion = m_scatt_pattern.memory_consumption ();
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      for (unsigned int h = 0; h < mp_state.get_n_groups (); ++h)
      {
        if (m_non_zero[g][h])
        {
          mem_consumtion += m_scatt_matrix[g][h].memory_consumption ();
        }
      }
    }
    return mem_consumtion;
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::vmult_add_g_h (BlockVector<double> & source,
                                             const BlockVector<double> & phi,
                                             const unsigned int g,
                                             const unsigned int h)
  {
    if (m_non_zero[g][h])
    {
      vmult_add_g_h_matrix_built (source, phi, g, h);
    }
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::vmult_add_g_h_matrix_built (BlockVector<double> & source,
                                                          const BlockVector<
                                                              double> & phi,
                                                          const unsigned int g,
                                                          const unsigned int h)
  {
    // the dofs of the scalar flux are at the beginning, so using the size of
    // the matrix to perform the product we are ignoring the dofs of the
    // boundaries
    for (unsigned int n = 0; n < m_scatt_matrix[g][h].n_block_rows (); ++n)
    {
      for (unsigned int m = 0; m < m_scatt_matrix[g][h].n_block_cols (); ++m)
      {
        m_scatt_matrix[g][h].block (n, m).vmult_add (source.block (n),
            phi.block (m));
      }
    }
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::setup_system ()
  {
    const unsigned int n_g = mp_state.get_n_groups ();
    const unsigned int n_m = mp_state.get_n_leg_mom ();
    m_scatt_pattern.reinit (n_m, n_m); /* Setting up scatt_pattern. */
    for (unsigned int i_ord = 0; i_ord < n_m; ++i_ord)
    {
      for (unsigned int j_ord = 0; j_ord < n_m; ++j_ord)
      {
        DynamicSparsityPattern csp;
        csp.reinit (mp_dof_handler.n_dofs (), mp_dof_handler.n_dofs ());
        DoFTools::make_flux_sparsity_pattern (mp_dof_handler, csp);
        csp.compress ();
        m_scatt_pattern.block (i_ord, j_ord).copy_from (csp);
      }
    }
    m_scatt_pattern.collect_sizes ();

    m_scatt_matrix.resize (n_g, std::vector<BlockSparseMatrix<double> > (n_g));
    for (unsigned int g = 0; g < n_g; ++g)
    {
      for (unsigned int h = 0; h < n_g; ++h)
      {
        if (m_non_zero[g][h])
        {
          m_scatt_matrix[g][h].reinit (m_scatt_pattern);
        }
      }
    }
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::build_system ()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;

    const unsigned int n_gauss_points = mp_dof_handler.get_fe ().degree + 1;
    info_box.initialize_gauss_quadrature (n_gauss_points, n_gauss_points,
        n_gauss_points);

    info_box.initialize_update_flags ();
    UpdateFlags update_flags = update_quadrature_points | update_values
                               | update_JxW_values; //update_JxW_values needed?
    info_box.add_update_flags (update_flags, true, true, true, true);

    info_box.initialize (mp_fe, mp_mapping);

    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    MeshWorker::DoFInfo<dim> dof_info_scalar (mp_dof_handler);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        if (m_non_zero[g][h])
        {
          for (unsigned int i_ord = 0; i_ord < n_leg_mom; ++i_ord)
          {
            for (unsigned int j_ord = 0; j_ord < n_leg_mom; ++j_ord)
            {
              MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
              assembler.initialize (m_scatt_matrix[g][h].block (i_ord, j_ord));

              MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
                  MeshWorker::IntegrationInfoBox<dim> > (
                  mp_dof_handler.begin_active (), mp_dof_handler.end (),
                  dof_info_scalar, info_box,
                  std_cxx11::bind (
                      &ScatteringMatrixBuilt<dim>::scatt_integrate_cell_term,
                      this, std_cxx11::_1, std_cxx11::_2, g, h), nullptr,
                  nullptr, assembler);
            }
          }
        }
      }
    }
  }

  template <int dim>
  void
  ScatteringMatrixBuilt<dim>::scatt_integrate_cell_term (DoFInfo &dinfo,
                                                         CellInfo &info,
                                                         const unsigned int g,
                                                         const unsigned int h)
  {
    FullMatrix<double> &local_matrix = dinfo.matrix (0).matrix;
    const FEValuesBase<dim> &fe_v = info.fe_values ();
    const std::vector<double> &JxW = fe_v.get_JxW_values ();

    const unsigned int mat_id = dinfo.cell->material_id ();
    /*
     * Here the sigmas should come from TH
     */
    double sigmas = mp_state.mp_data.mp_mat.m_xs.at (mat_id).m_sigmas[g][h];

    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          local_matrix (i, j) += sigmas * fe_v.shape_value (i, q) * JxW[q]
                                 * fe_v.shape_value (j, q);
        }
      }
    }
  }

  template class ScatteringMatrixBuilt<1> ;
  template class ScatteringMatrixBuilt<2> ;
  template class ScatteringMatrixBuilt<3> ;

} // end of namespace Forest

