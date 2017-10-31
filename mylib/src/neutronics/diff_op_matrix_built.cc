/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op_matrix_built.cc
 * @brief  Implementation of class template DiffOpMatrixBuilt
 */

#include "neutronics/diff_op_matrix_built.h"

#include "utils/forest_utils_dealii.h"
#include "geometry/path_to_gridgenerator_dealii.h"
#include "geometry/prob_geom.h"
#include "input/input.h"
#include "neutronics/state.h"

#include "deal.II/base/config.h"         // for DEAL_II_VERSION_GTE
#include "deal.II/base/geometry_info.h"  // for GeometryInfo
#include "deal.II/base/point.h"          // for Point
#include "deal.II/base/tensor.h"         // for Tensor
#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/dofs/dof_tools.h"      // for FiniteElement, Mapping
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc
#include "deal.II/lac/vector.h"          // for Vector
#include "deal.II/lac/full_matrix.h"     // for FullMatrix
#include "deal.II/fe/fe_values.h"        // for FEFaceValues, etc

#include "deal.II/meshworker/dof_info.h"
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/assembler.h"
#include "deal.II/meshworker/loop.h"
#include "deal.II/integrators/laplace.h"
#include "deal.II/lac/sparse_matrix.h"

#include "deal.II/lac/precondition.h"
#include "deal.II/lac/precondition_block.h"
#include "deal.II/lac/sparse_mic.h"
#include "deal.II/lac/sparse_ilu.h"

#include <fstream>
#include <utility>                       // for pair
#include <vector>                        // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  DiffOpMatrixBuilt<dim>::DiffOpMatrixBuilt (State<dim> & state)
      : mp_state (state),
        mp_fe (state.get_fe ()),
        mp_mapping (state.get_mapping ()),
        mp_dof_handler (mp_state.get_dof_handler ())
  {
    this->setup_system ();
    this->assemble_system ();

    this->setup_prec ();
    this->assemble_prec ();
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::vmult (BlockVector<double> & phi,
                                 const BlockVector<double> & source,
                                 const unsigned int g)
  {
    this->vmult_matrix_built_g_h (phi, source, g);
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::prec (BlockVector<double> & phi,
                                const BlockVector<double> & source,
                                const unsigned int g) const
  {
    m_prec[g].vmult (phi.block (0), source.block (0));
    /*
    bool print_matrix = false;
    if (print_matrix)
    {
      std::ofstream out ("matrix.dat");
      const unsigned int precision = 6;
      const bool scientific = true;
      const unsigned int width = 6;
      const char *zero_string = "0.0000";
      const double denominator = 1.;
      m_Diff_matrix[g].block (0, 0).print_formatted (out, precision, scientific,
          width, zero_string, denominator);
      Assert(false, ExcMessage("Stop here!"));
    }

    const unsigned int prec = 2;
    if (prec == 0)
    {
      SparseMIC<double> precondition;
      const double strengthen_diagonal=1;
      const unsigned int extra_off_diagonals=0;
      precondition.initialize (m_Diff_matrix[g].block (0, 0),
          SparseMIC<double>::AdditionalData(strengthen_diagonal,extra_off_diagonals));
      precondition.vmult (phi.block (0), source.block (0));
    }
    else if (prec == 1)
    {
      SparseILU<double> precondition;
      precondition.initialize (m_Diff_matrix[g].block (0, 0),
        SparseILU<double>::AdditionalData());
      precondition.vmult (phi.block (0), source.block (0));
    }
    else if (prec == 1)
    {
      PreconditionBlockSSOR<SparseMatrix<double>, double> precondition;
      const unsigned int dofs_per_cell = mp_state.get_fe ().dofs_per_cell;
      precondition.initialize (m_Diff_matrix[g].block (0, 0), dofs_per_cell);
      precondition.vmult (phi.block (0), source.block (0));
    }
    else if (prec == 2)
    {
      PreconditionBlockJacobi<SparseMatrix<double>, double> precondition;
      const unsigned int dofs_per_cell = mp_state.get_fe ().dofs_per_cell;
      precondition.initialize (m_Diff_matrix[g].block (0, 0), dofs_per_cell);
      precondition.vmult (phi.block (0), source.block (0));
    }
    else if (prec == 3)
    {
      PreconditionSSOR<SparseMatrix<double> > precondition;
      precondition.initialize (m_Diff_matrix[g].block (0, 0), 1.);
      precondition.vmult (phi.block (0), source.block (0));
    }
    else if (prec == 4)
    {
      PreconditionJacobi<SparseMatrix<double> > precondition;
      precondition.initialize (m_Diff_matrix[g].block (0, 0), .6);
      precondition.vmult (phi.block (0), source.block (0));

    }
    else if (prec == 5)
    {
      //const unsigned int n_m = 1;
      for (unsigned int i = 0; i < mp_dof_handler.n_dofs (); ++i)
        phi.block (0)[i] = source.block (0)[i] * m_prec_matrix[g].block (0)[i];
    }
    else
    {
      phi = source;
    }*/

  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::vmult_matrix_built_g_h (BlockVector<double> & dst,
                                                  const BlockVector<double> & src,
                                                  const unsigned int g)
  {
    /*std::cout << "DiffOpMatrixBuilt<dim>::vmult_matrix_free_g_h" << std::endl;
    std::cout << "dst.n_blocks() = " << dst.n_blocks() << std::endl;
    std::cout << "src.n_blocks() = " << src.n_blocks() << std::endl;*/
    //m_Diff_matrix[g].vmult_add (dst, src);
    //dst = 0.;
    const unsigned int n_m = 1;
    for (unsigned int m = 0; m < n_m; ++m)
    {
      for (unsigned int l = 0; l < n_m; ++l)
      {
        m_Diff_matrix[g].block (m, l).vmult_add (dst.block (m), src.block (l));
      }
    }
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::setup_system ()
  {
    const unsigned int n_g = mp_state.get_n_groups ();
    const unsigned int n_m = 1;
    m_Diff_pattern.reinit (n_m, n_m); /* Setting up scatt_pattern. */
    for (unsigned int m = 0; m < n_m; ++m)
    {
      for (unsigned int l = 0; l < n_m; ++l)
      {
        DynamicSparsityPattern csp;
        csp.reinit (mp_dof_handler.n_dofs (), mp_dof_handler.n_dofs ());
        DoFTools::make_flux_sparsity_pattern (mp_dof_handler, csp);
        csp.compress ();
        m_Diff_pattern.block (m, l).copy_from (csp);
      }
    }
    m_Diff_pattern.collect_sizes ();

    m_Diff_matrix.resize (n_g);
    for (unsigned int g = 0; g < n_g; ++g)
    {
      m_Diff_matrix[g].reinit (m_Diff_pattern);
    }
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::assemble_system ()
  {

    MeshWorker::IntegrationInfoBox<dim> info_box;

    const unsigned int n_gauss_points = mp_dof_handler.get_fe ().degree + 1;
    info_box.initialize_gauss_quadrature (n_gauss_points, n_gauss_points,
        n_gauss_points);

    info_box.initialize_update_flags ();
    UpdateFlags update_flags = update_quadrature_points | update_values
                               | update_gradients | update_JxW_values; //update_JxW_values needed?
    info_box.add_update_flags (update_flags, true, true, true, true);
    info_box.initialize (mp_fe, mp_mapping);

    const unsigned int n_groups = mp_state.get_n_groups ();
    // pure diffusion with n_mom = 1 (for higher order we get SPn)
    const unsigned int n_mom = 1;

    MeshWorker::DoFInfo<dim> dof_info_scalar (mp_dof_handler);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      bool only_mass = false;
      bool use_cell = true;
      bool use_boundary = true;
      bool use_face = true;
      MatrixIntegrator<dim> integrator_Diff (use_cell, use_boundary, use_face,
          g, g, mp_state.mp_data, only_mass);
      for (unsigned int m = 0; m < n_mom; ++m)
      {
        for (unsigned int l = 0; l < n_mom; ++l)
        {
          MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
          assembler.initialize (m_Diff_matrix[g].block (m, l));
          MeshWorker::integration_loop<dim, dim> (
              mp_dof_handler.begin_active (), mp_dof_handler.end (),
              dof_info_scalar, info_box, integrator_Diff, assembler);
        }
      }
      m_Diff_matrix[g].compress (VectorOperation::add);
    }
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::setup_prec ()
  {
    const unsigned int n_g = mp_state.get_n_groups ();
    /*const unsigned int n_m = 1;*/
    m_prec.resize (n_g);
  }

  template <int dim>
  void
  DiffOpMatrixBuilt<dim>::assemble_prec ()
  {
    const unsigned int n_g = mp_state.get_n_groups ();
    //const unsigned int n_m = 1;
    for (unsigned int g = 0; g < n_g; ++g)
    {
      //m_prec[g].initialize (m_Diff_matrix[g].block (0, 0), 1.0);
      m_prec[g].initialize (m_Diff_matrix[g].block (0, 0),
                            mp_state.get_fe ().dofs_per_cell);
      //m_prec[g].initialize (m_Diff_matrix[g].block (0, 0), SparseILU<double>::AdditionalData());
      //m_prec[g].initialize (m_Diff_matrix[g].block (0, 0), SparseMIC<double>::AdditionalData());
    }
  }

  template <int dim>
  unsigned int
  DiffOpMatrixBuilt<dim>::get_n_groups () const
  {
    return mp_state.get_n_groups ();
  }

  template <int dim>
  double
  DiffOpMatrixBuilt<dim>::memory_consumption () const
  {
    const unsigned int n_groups = mp_state.get_n_groups ();
    double memory_consumption = 0;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      memory_consumption += m_Diff_matrix[g].memory_consumption ();
    }
    return memory_consumption;
  }

  template class DiffOpMatrixBuilt<1> ;
  template class DiffOpMatrixBuilt<2> ;
  template class DiffOpMatrixBuilt<3> ;

} // end of namespace Forest
