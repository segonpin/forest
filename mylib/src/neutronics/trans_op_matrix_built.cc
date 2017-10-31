/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op_matrix_built.cc
 * @brief  Implementation of class template TransOpMatrixBuilt
 */

#include "neutronics/trans_op_matrix_built.h"  // for TransOpMatrixBuilt, etc

#include "neutronics/state.h"
#include "input/input.h"
#include "input/input_geom.h"
#include "input/input_mat.h"
#include "angle/quadrature_base.h"
#include "utils/forest_utils_dealii.h" // for vector_to_point
#include "neutronics/downstream_new.h"             // for compute_downstream_new
#include "neutronics/path_to_intergridmap_dealii.h"
#include "geometry/prob_geom.h"

#include "deal.II/base/table.h"         // for Table
#include "deal.II/base/std_cxx11/bind.h"
#include "deal.II/base/point.h"         // for Point
#include "deal.II/base/tensor.h"        // for Tensor
#include "deal.II/base/types.h"         // for global_dof_index

#include "deal.II/lac/full_matrix.h"    // for FullMatrix
#include "deal.II/lac/sparse_matrix.h"  // for SparseMatrix
#include "deal.II/lac/constraint_matrix.h"  // for ConstraintMatrix
#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/precondition_block.h"  // for PreconditionBlockSOR
#include "deal.II/lac/vector.h"         // for Vector

#include "deal.II/dofs/dof_tools.h"     // for make_flux_sparsity_pattern

#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc
#include "deal.II/fe/fe_values.h"       // for FEValuesBase

#include "deal.II/meshworker/dof_info.h"  // for DoFInfo
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/loop.h"    // for loop
#include "deal.II/meshworker/simple.h"  // for MatrixSimple

#include <vector>                       // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  TransOpMatrixBuilt<dim>::TransOpMatrixBuilt (State<dim> & state,
                                               std::shared_ptr<
                                                   QuadratureBase<dim> > & ord)
      : mp_state (state),
        mp_ord (ord),
        mp_dof_handler (mp_state.get_dof_handler ()),
        m_dof_handler_angular (get_n_angles ())
  {
    m_projected_block_dst.resize (get_n_angles ());
    m_projected_block_src.resize (get_n_angles ());
    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      m_projected_block_dst[i_ord].reinit (mp_dof_handler.n_dofs ());
      m_projected_block_src[i_ord].reinit (mp_dof_handler.n_dofs ());
    }
    build ();
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::build ()
  {
    setup_angular_dof_handler ();
    setup_system ();
    build_trans_matrix ();
    build_prec_matrix ();
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::apply_L0_inv_g_iord (Vector<double> & phi,
                                                const Vector<double> & source,
                                                const unsigned int g,
                                                const unsigned int i_ord)
  {
    apply_L0_inv_matrix_built (phi, source, g, i_ord);
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::setup_system ()
  {
    m_trans_pattern.resize (get_n_angles ());
    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      DynamicSparsityPattern csp;
      csp.reinit (m_dof_handler_angular[i_ord].n_dofs (),
          m_dof_handler_angular[i_ord].n_dofs ());
      DoFTools::make_flux_sparsity_pattern (m_dof_handler_angular[i_ord], csp);
      m_trans_pattern[i_ord].copy_from (csp);
    }

    m_trans_matrix.resize (get_n_groups (),
        std::vector<SparseMatrix<double> > (get_n_angles ()));

    m_prec_matrix.resize (get_n_groups (),
        std::vector<PreconditionBlock<SparseMatrix<double>, double> > (
            get_n_angles ()));
    /*m_prec_matrix.resize (get_n_groups (),
     std::vector<PreconditionBlockSOR<SparseMatrix<double>, double> > (
     get_n_angles ()));*/
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::setup_angular_dof_handler ()
  {
    m_renumbering_dof_indices.resize (get_n_angles (),
        std::vector<types::global_dof_index> (mp_dof_handler.n_dofs ()));
    m_reverse_dof_indices.resize (get_n_angles (),
        std::vector<types::global_dof_index> (mp_dof_handler.n_dofs ()));

    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      m_dof_handler_angular[i_ord].initialize (mp_state.mp_geom.m_tria,
          mp_state.get_fe ());
      const Point<dim> down_direction = vector_to_point<dim> (
          mp_ord->get_q (i_ord));
      DoFRenumbering::compute_downstream_new (m_renumbering_dof_indices[i_ord],
          m_reverse_dof_indices[i_ord], m_dof_handler_angular[i_ord],
          down_direction);
      m_dof_handler_angular[i_ord].renumber_dofs (
          m_renumbering_dof_indices[i_ord]);
    }

    // resize the containers
    m_grid_m_to_d_map.resize (get_n_angles ());
    m_grid_d_to_m_map.resize (get_n_angles ());
    // We also initialize the mapping between the different orderings.
    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      m_grid_m_to_d_map[i_ord].make_mapping (mp_dof_handler,
          m_dof_handler_angular[i_ord]);
      m_grid_d_to_m_map[i_ord].make_mapping (m_dof_handler_angular[i_ord],
          mp_dof_handler);
    }

    // We need the constraints for the grid_map, even if we keep them
    // uninitialized
    constraints.clear ();
    constraints.close ();
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::build_trans_matrix ()
  {
    for (unsigned int g = 0; g < get_n_groups (); ++g)
    {
      for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
      {
        m_trans_matrix[g][i_ord].reinit (m_trans_pattern[i_ord]);
      }
    }

    MeshWorker::IntegrationInfoBox<dim> info_box;

    const unsigned int n_gauss_points = mp_state.get_fe ().degree + 1;
    info_box.initialize_gauss_quadrature (n_gauss_points, n_gauss_points,
        n_gauss_points);

    info_box.initialize_update_flags ();
    UpdateFlags update_flags = update_quadrature_points | update_values
                               | update_gradients | update_JxW_values;
    info_box.add_update_flags (update_flags, true, true, true, true);
    info_box.initialize (mp_state.get_fe (), mp_state.get_mapping ());

    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      // updating the streaming direction
      const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));
      const double weight = mp_ord->get_w (i_ord);
      for (unsigned int g = 0; g < get_n_groups (); ++g)
      {
        MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
        assembler.initialize (m_trans_matrix[g][i_ord]);
        MeshWorker::DoFInfo<dim> dof_info (m_dof_handler_angular[i_ord]);
        MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
            MeshWorker::IntegrationInfoBox<dim> > (
            m_dof_handler_angular[i_ord].begin_active (),
            m_dof_handler_angular[i_ord].end (), dof_info, info_box,
            std_cxx11::bind (
                &TransOpMatrixBuilt<dim>::trans_integrate_cell_term, this,
                std_cxx11::_1, std_cxx11::_2, g, beta, weight),
            std_cxx11::bind (
                &TransOpMatrixBuilt<dim>::integrate_outflux_boundary_term, this,
                std_cxx11::_1, std_cxx11::_2, beta),
            std_cxx11::bind (&TransOpMatrixBuilt<dim>::integrate_face_term,
                this, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3,
                std_cxx11::_4, beta), assembler);
      } // g
    } // i_ord
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::build_prec_matrix ()
  {
    for (unsigned int g = 0; g < get_n_groups (); ++g)
    {
      for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
      {
        m_prec_matrix[g][i_ord].initialize (m_trans_matrix[g][i_ord],
            mp_state.get_fe ().dofs_per_cell);
        //m_prec_matrix[g][i_ord].invert_diagblocks ();
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::apply_L0_inv_matrix_built (Vector<double> & flux,
                                                      const Vector<double> & source,
                                                      const unsigned int g,
                                                      const unsigned int i_ord)
  {
    // first we project the mesh from m_dof_handler_angular to mp_dof_handler,
    // and then we solve the system, where the system is built with the order
    // in m_dof_handler_angular std::cout << "project_scalar_to_angular_mesh\n";
    project_scalar_to_angular_mesh (i_ord, source,
        m_projected_block_src[i_ord]);
    //deallog << "sweep for (g,i_ord) = (" << g << "," << i_ord << ")" << std::endl;
    // Then we solve
    //m_prec_matrix[g][i_ord].backward_step(m_projected_block_dst[i_ord],m_projected_block_dst[i_ord], m_projected_block_src[i_ord],false);
    m_prec_matrix[g][i_ord].forward_step (m_projected_block_dst[i_ord],
        m_projected_block_dst[i_ord], m_projected_block_src[i_ord], false);
    //m_prec_matrix[g][i_ord].vmult(m_projected_block_dst[i_ord], m_projected_block_src[i_ord]);
    /*m_solver.solve (m_trans_matrix[g][i_ord], m_projected_block_dst,
     m_projected_block_src, m_prec_matrix[g][i_ord]);
     deallog << "n_iter = " << m_solver_control.last_step() << std::endl;*/
    // then we project back
    project_angular_to_scalar_mesh (i_ord, m_projected_block_dst[i_ord], flux);
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::trans_integrate_cell_term (DoFInfo &dinfo,
                                                      CellInfo &info,
                                                      const unsigned int g,
                                                      const Point<dim> & beta,
                                                      const double /* weight */)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values ();
    FullMatrix<double> &local_matrix = dinfo.matrix (0).matrix;
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    /*
     * Here the sigmat should come from TH
     */
    const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
        dinfo.cell->material_id ()).m_sigmat[g];
    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          local_matrix (i, j) += -beta * fe_v.shape_grad (i, q)
                                 * fe_v.shape_value (j, q) * JxW[q]
                                 + sigmat * fe_v.shape_value (i, q)
                                   * fe_v.shape_value (j, q) * JxW[q];
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::integrate_outflux_boundary_term (DoFInfo &dinfo,
                                                            CellInfo &info,
                                                            const Point<dim> &beta)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values ();
    FullMatrix<double> &local_matrix = dinfo.matrix (0).matrix;

    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Tensor<1, dim> > &normals =
        fe_v.get_all_normal_vectors ();

    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      const double beta_n = beta * normals[q];
      if (beta_n > 0)
      {
        for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
          {
            local_matrix (i, j) += beta_n * fe_v.shape_value (j, q)
                                   * fe_v.shape_value (i, q) * JxW[q];
          }
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::integrate_influx_boundary_term (DoFInfo &dinfo,
                                                           CellInfo &info,
                                                           const Point<dim> & beta,
                                                           const unsigned int boundary_face)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values ();
    FullMatrix<double> &local_matrix = dinfo.matrix (0).matrix;
    //Vector<double> &local_vector = dinfo.vector(0).block(0);

    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Tensor<1, dim> > &normals =
        fe_v.get_all_normal_vectors ();
    const unsigned int bid = dinfo.face->boundary_id ();

    if (mp_state.mp_data.mp_geom.m_core.m_bcs[bid] == 2 and bid == boundary_face)
    {
      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
      {
        const double beta_n = beta * normals[q];
        if (beta_n < 0)
        {
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
            {
              local_matrix (i, j) += beta_n * fe_v.shape_value (j, q)
                                     * fe_v.shape_value (i, q) * JxW[q];
            }
          }
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::integrate_face_term (DoFInfo &dinfo1,
                                                DoFInfo &dinfo2,
                                                CellInfo &info1,
                                                CellInfo &info2,
                                                const Point<dim> &beta)
  {
    const FEValuesBase<dim> &fe_v = info1.fe_values ();
    const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values ();

    FullMatrix<double> &u1_v1_matrix = dinfo1.matrix (0, false).matrix;
    FullMatrix<double> &u2_v1_matrix = dinfo1.matrix (0, true).matrix;
    FullMatrix<double> &u1_v2_matrix = dinfo2.matrix (0, true).matrix;
    FullMatrix<double> &u2_v2_matrix = dinfo2.matrix (0, false).matrix;

    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    const std::vector<Tensor<1, dim> > &normals =
        fe_v.get_all_normal_vectors ();

    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      const double beta_n = beta * normals[q];
      if (beta_n > 0)
      {
        for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
          {
            u1_v1_matrix (i, j) += beta_n * fe_v.shape_value (j, q)
                                   * fe_v.shape_value (i, q) * JxW[q];
          }
        }
        for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        {
          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
          {
            u1_v2_matrix (k, j) -= beta_n * fe_v.shape_value (j, q)
                                   * fe_v_neighbor.shape_value (k, q) * JxW[q];
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        {
          for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
          {
            u2_v1_matrix (i, l) += beta_n * fe_v_neighbor.shape_value (l, q)
                                   * fe_v.shape_value (i, q) * JxW[q];
          }
        }

        for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        {
          for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
          {
            u2_v2_matrix (k, l) -= beta_n * fe_v_neighbor.shape_value (l, q)
                                   * fe_v_neighbor.shape_value (k, q) * JxW[q];
          }
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::project_scalar_to_angular_mesh (const unsigned int i,
                                                           const Vector<double> & vec_scalar,
                                                           Vector<double> & vec_angular)
  {
    VectorTools::interpolate_to_different_mesh (m_grid_m_to_d_map[i],
        vec_scalar, constraints, vec_angular);
  }

  template <int dim>
  void
  TransOpMatrixBuilt<dim>::project_angular_to_scalar_mesh (const unsigned int i,
                                                           const Vector<double> & vec_angular,
                                                           Vector<double> & vec_scalar)
  {
    VectorTools::interpolate_to_different_mesh (m_grid_d_to_m_map[i],
        vec_angular, constraints, vec_scalar);
  }

  template <int dim>
  unsigned int
  TransOpMatrixBuilt<dim>::get_n_groups () const
  {
    return mp_state.get_n_groups ();
  }

  template <int dim>
  unsigned int
  TransOpMatrixBuilt<dim>::get_n_angles () const
  {
    return mp_ord->get_n_angles ();
  }

  template <int dim>
  double
  TransOpMatrixBuilt<dim>::memory_consumption () const
  {
    double memory_consumption = 0;

    memory_consumption += (m_projected_block_dst[0].memory_consumption ()
        + m_projected_block_dst[0].memory_consumption ())
                          * get_n_angles ();

    for (unsigned int g = 0; g < get_n_groups (); ++g)
      for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
        memory_consumption += m_trans_matrix[g][i_ord].memory_consumption ();

    for (unsigned int g = 0; g < get_n_groups (); ++g)
      for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
        memory_consumption += m_prec_matrix[g][i_ord].memory_consumption ();

    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
      memory_consumption += m_dof_handler_angular[i_ord].memory_consumption ();

    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
      memory_consumption += m_grid_m_to_d_map[i_ord].memory_consumption ();

    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
      memory_consumption += m_grid_d_to_m_map[i_ord].memory_consumption ();

    return memory_consumption;
  }

  template class TransOpMatrixBuilt<1> ;
  template class TransOpMatrixBuilt<2> ;
  template class TransOpMatrixBuilt<3> ;

} // end of namespace Forest
