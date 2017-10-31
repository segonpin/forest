/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_matrix_built.cc
 * @brief  Implementation of class template FissionSourceMatrixBuilt
 */

#include "neutronics/fission_source_matrix_built.h"

#include "neutronics/state.h"
#include "input/input.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc
#include "deal.II/fe/mapping_q1.h"       // for MappingQ1
#include "deal.II/lac/dynamic_sparsity_pattern.h"  // for DynamicSparsityPat...
#include "deal.II/lac/full_matrix.h"
#include "deal.II/lac/sparse_matrix.h"

#include "deal.II/meshworker/dof_info.h"
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/simple.h"
#include "deal.II/meshworker/loop.h"
#include "deal.II/meshworker/local_integrator.h"   // for DoFInfo
#include "deal.II/meshworker/vector_selector.h"    // for FEValuesBase
#include "deal.II/base/std_cxx11/bind.h"

namespace Forest
{
  using namespace dealii;

  template <int dim>
  FissionSourceMatrixBuilt<dim>::FissionSourceMatrixBuilt (State<dim> & state)
      : mp_state (state),
        mp_dof_handler (state.get_dof_handler ())
  {
    // We first setup the non zero pattern, in order to avoid operations later
    this->set_non_zero_pattern ();
    // Now we build the matrix
    this->build ();
  }

  template <int dim>
  void
  FissionSourceMatrixBuilt<dim>::vmult (SnVector & source,
                                        const SnVector & phi)
  {
    source.m_data[0].block (0) = 0.0;
    vmult_add (source, phi);
  }

  template <int dim>
  void
  FissionSourceMatrixBuilt<dim>::vmult_add (SnVector & source,
                                            const SnVector & phi)
  {
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      if (m_non_zero[g])
      {
        m_nusigf_matrix[g].vmult_add (source.m_data[0].block (0),
            phi.m_data[g].block (0));
      }
    }
  }

  template <int dim>
  void
  FissionSourceMatrixBuilt<dim>::set_non_zero_pattern ()
  {
    // we consider zero anything below this value;
    const double eps = 1.e-14;
    const unsigned int n_groups = mp_state.get_n_groups ();
    m_non_zero.resize (n_groups, false);
    double nusigf;
    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      for (cell = mp_dof_handler.begin_active (); cell != endc; ++cell)
      {
        /*
         * Here the nusigf should come from TH
         */
        nusigf =
            mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_nusigf[g];
        if (std::abs (nusigf) > eps)
        {
          m_non_zero[g] = true;
          // Once we find a non-zero value, then stop the loop for this group
          continue;
        }
      }
    }
  }

  template <int dim>
  void
  FissionSourceMatrixBuilt<dim>::build ()
  {
    DynamicSparsityPattern csp;
    csp.reinit (mp_dof_handler.n_dofs (), mp_dof_handler.n_dofs ());
    DoFTools::make_flux_sparsity_pattern (mp_dof_handler, csp);
    csp.compress ();

    m_nusigf_pattern.copy_from (csp);

    m_nusigf_matrix.resize (mp_state.get_n_groups ());
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      if (m_non_zero[g])
      {
        m_nusigf_matrix[g].reinit (m_nusigf_pattern);
      }
    }

    MeshWorker::IntegrationInfoBox<dim> info_box;

    const unsigned int n_gauss_points = mp_dof_handler.get_fe ().degree + 1;
    info_box.initialize_gauss_quadrature (n_gauss_points, n_gauss_points,
        n_gauss_points);
    info_box.initialize_update_flags ();
    UpdateFlags update_flags = update_quadrature_points | update_values
                               | update_JxW_values;
    info_box.add_update_flags (update_flags, true, true, true, true);

    const FiniteElement<dim> & fe (mp_state.get_fe ());
    const MappingQ1<dim> & mapping (mp_state.get_mapping ());
    info_box.initialize (fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info_scalar (mp_dof_handler);
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      if (m_non_zero[g])
      {
        MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
        assembler.initialize (m_nusigf_matrix[g]);
        MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
            MeshWorker::IntegrationInfoBox<dim> > (
            mp_dof_handler.begin_active (), mp_dof_handler.end (),
            dof_info_scalar, info_box,
            std_cxx11::bind (
                &FissionSourceMatrixBuilt<dim>::nsigf_integrate_cell_term, this,
                std_cxx11::_1, std_cxx11::_2, g), nullptr, nullptr, assembler);
      }
    }
  }

  template <int dim>
  void
  FissionSourceMatrixBuilt<dim>::nsigf_integrate_cell_term (MeshWorker::DoFInfo<
                                                                dim> &dinfo,
                                                            MeshWorker::IntegrationInfo<
                                                                dim> &info,
                                                            const unsigned int g)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values ();
    FullMatrix<double> &local_matrix = dinfo.matrix (0).matrix;
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    /*
     * Here the nusigf should come from TH
     */
    double nusigf =
        mp_state.mp_data.mp_mat.m_xs.at (dinfo.cell->material_id ()).m_nusigf[g];

    for (unsigned int q_point = 0; q_point < fe_v.n_quadrature_points;
        ++q_point)
    {
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          local_matrix (i, j) +=
              (nusigf * fe_v.shape_value (i, q_point)
               * fe_v.shape_value (j, q_point) * JxW[q_point]);
        }
      }
    }
  }

  template <int dim>
  double
  FissionSourceMatrixBuilt<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption += m_nusigf_pattern.memory_consumption ();
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      if (m_non_zero[g])
      {
        memory_consumption += m_nusigf_matrix[g].memory_consumption ();
      }
    }
    return memory_consumption;
  }

  template class FissionSourceMatrixBuilt<1> ;
  template class FissionSourceMatrixBuilt<2> ;
  template class FissionSourceMatrixBuilt<3> ;

} // end of namespace Forest
