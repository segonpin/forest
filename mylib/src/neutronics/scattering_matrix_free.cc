/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_matrix_free.cc
 * @brief  Implementation of class template ScatteringMatrixFree
 */

#include "neutronics/scattering_matrix_free.h"

#include "neutronics/state.h"
#include "geometry/prob_geom.h"
#include "input/input.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/base/quadrature_lib.h"  // for QGauss
#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/fe/fe.h"
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc
#include "deal.II/fe/fe_values.h"
#include "deal.II/fe/mapping.h"
#include "deal.II/grid/tria.h"          // for Triangulation
#include "deal.II/lac/block_vector_base.h"
#include "deal.II/lac/vector.h"         // for Vector
#include "deal.II/grid/grid_generator.h"

#include <vector>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  ScatteringMatrixFree<dim>::ScatteringMatrixFree (State<dim> & state)
      : mp_state (state),
        mp_fe (state.get_fe ()),
        mp_mapping (state.get_mapping ()),
        mp_dof_handler (state.get_dof_handler ())
  {
    // We first setup the non zero pattern, in order to avoid operations later
    this->set_non_zero_pattern ();
    //this->print_non_zero_pattern();
  }

  template <int dim>
  void
  ScatteringMatrixFree<dim>::set_non_zero_pattern ()
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
  ScatteringMatrixFree<dim>::print_non_zero_pattern ()
  {
    for (unsigned int g = 0; g < m_non_zero.size(); ++g)
    {
      for (unsigned int h = 0; h < m_non_zero[g].size(); ++h)
      {
        dealii::deallog << " " << m_non_zero[g][h];
      }
      dealii::deallog << std::endl;
    }
    dealii::deallog << std::endl;
    //sleep(0.01);
  }

  template <int dim>
  bool
  ScatteringMatrixFree<dim>::no_upscatt(const unsigned int g) const
  {
     for (unsigned int h = g+1; h < m_non_zero[g].size(); ++h)
     {
       if (m_non_zero[g][h])
       {
         return false;
       }
     }
     return true;
  }

  template <int dim>
  void
  ScatteringMatrixFree<dim>::vmult_add_g_h (BlockVector<double> & source,
                                            const BlockVector<double> & phi,
                                            const unsigned int g,
                                            const unsigned int h)
  {
    if (m_non_zero[g][h])
    {
      vmult_add_g_h_matrix_free(source, phi, g, h);
    }
  }

  template <int dim>
  double
  ScatteringMatrixFree<dim>::memory_consumption () const
  {
    double mem_consumtion = 0;
    return mem_consumtion;
  }


  template <int dim>
  void
  ScatteringMatrixFree<dim>::vmult_add_g_h_matrix_free (BlockVector<double> & source,
                                            const BlockVector<double> & phi,
                                            const unsigned int g,
                                            const unsigned int h)
  {
    const FiniteElement<dim> & fe (mp_state.get_fe ());
    const Mapping<dim> & mapping (mp_state.get_mapping ());
    QGauss<dim> quadrature_formula (fe.degree + 1);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size ();
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell);
    Vector<double> tmp_vec (n_q_points);

    FEValues<dim> fe_values_reference (mapping, fe, quadrature_formula,
        update_values);
    Triangulation<dim, dim> reference_cell;
    GridGenerator::hyper_cube (reference_cell, 0., 1.);
    fe_values_reference.reinit (reference_cell.begin ());
    FEValues<dim> fe_values (fe, quadrature_formula,
        update_JxW_values);

    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      /*
       * Here the sigmas should come from TH
       */
      double sigmas =
          mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_sigmas[g][h];

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        cell_src (i) = phi (local_dof_indices[i]);
      }

      tmp_vec = 0;
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          tmp_vec (q) += fe_values_reference.shape_value (i, q) * cell_src (i);
        }
      }

      // multiply by coefficient and integration weight
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        tmp_vec (q) *= fe_values.JxW (q) * sigmas;
      }

      cell_dst = 0;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          cell_dst (i) += fe_values_reference.shape_value (i, q) * tmp_vec (q);
        }
      }

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        source (local_dof_indices[i]) += cell_dst (i);
      }
    }
  }

//  // This is not working yet, and I do not know why
//
//  template <int dim>
//  template <int fe_degree>
//  void
//  ScatteringMatrixFree<dim>::vmult_add_g_h_k (BlockVector<double> & source,
//                                              const BlockVector<double> & phi,
//                                              const unsigned int g,
//                                              const unsigned int h)
//  {
//    //const FiniteElement<dim> & fe (mp_state.get_fe ());
//
//    const unsigned int dofs_per_cell = mp_state.get_fe ().dofs_per_cell;
//    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
//    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell);
//    UpdateFlags flags = update_values | update_JxW_values;
//
//    std::cout << "bug 1" << std::endl;
//    FEEvaluation<dim, fe_degree> fe_eval (mp_state.get_fe (),
//        QGauss<1> (fe_degree + 1), flags);
//    std::cout << "bug 2" << std::endl;
//
//    double sigmas;
//    typename DoFHandler<dim>::active_cell_iterator cell =
//        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
//    for (; cell != endc; ++cell)
//    {
//      sigmas =
//          mp_state.mp_data.mp_mat.xs.at (cell->material_id ()).sigmas[g][h];
//      fe_eval.reinit (cell);
//
//      cell->get_dof_indices (local_dof_indices);
//      for (unsigned int i = 0; i < dofs_per_cell; ++i)
//        cell_src (i) = phi (local_dof_indices[i]);
//      fe_eval.read_dof_values (cell_src);
//
//      fe_eval.evaluate (true, false, false);
//      for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
//        fe_eval.submit_value (sigmas * fe_eval.get_value (q), q);
//
//      fe_eval.integrate (true, false);
//      fe_eval.distribute_local_to_global (cell_dst);
//
//      for (unsigned int i = 0; i < dofs_per_cell; ++i)
//        source (local_dof_indices[i]) += cell_dst (i);
//    }
//  }
//
//  template <int dim>
//  void
//  ScatteringMatrixFree<dim>::vmult_add_g_h (BlockVector<double> & source,
//                                            const BlockVector<double> & phi,
//                                            const unsigned int g,
//                                            const unsigned int h)
//  {
//    unsigned int fe_degree = mp_state.mp_data.mp_settings.get_fe_degree ();
//    switch (fe_degree)
//    {
//      case 1:
//        vmult_add_g_h_k<1> (source, phi, g, h);
//        break;
//      case 2:
//        vmult_add_g_h_k<2> (source, phi, g, h);
//        break;
//      case 3:
//        vmult_add_g_h_k<3> (source, phi, g, h);
//        break;
//      case 4:
//        vmult_add_g_h_k<4> (source, phi, g, h);
//        break;
//      case 5:
//        vmult_add_g_h_k<5> (source, phi, g, h);
//        break;
//      case 6:
//        vmult_add_g_h_k<6> (source, phi, g, h);
//        break;
//      case 7:
//        vmult_add_g_h_k<7> (source, phi, g, h);
//        break;
//      default:
//        Assert(false, ExcMessage("fe_degree not allowed."));
//        break;
//    }
//  }

//  //------------------------------------------------------------------------------
//  //  This should wait until next update of my version of the library
//  //------------------------------------------------------------------------------
//  template <int dim>
//  void
//  ScatteringMatrixFree<dim>::set_coeff (const unsigned int g,
//                                        const unsigned int h)
//  {
//    /* Resize the vector and initialize to zero. */
//    coeff.resize (mp_state.mp_geom.tria.n_active_cells (), 0.);
//    /* Add the values for the groups we are calculating. */
//    typename DoFHandler<dim, dim>::active_cell_iterator cell =
//        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
//    for (; cell != endc; ++cell)
//    {
//      const unsigned int user_index = cell->user_index ();
//      const unsigned int mat = cell->material_id ();
//      coeff[user_index] = mp_state.mp_data.mp_mat.xs.at (mat).sigmas[g][h];
//    }
//  }
//
//  template <int dim>
//  double
//  ScatteringMatrixFree<dim>::get_coeff (const unsigned int user_index)
//  {
//    return coeff[user_index];
//  }
//
//  template <int dim>
//  template <int fe_degree, typename number>
//  void
//  ScatteringMatrixFree<dim>::local_apply (const dealii::MatrixFree<dim, number> &data,
//                                          Vector<double> &dst,
//                                          const Vector<double> &src,
//                                          const std::pair<unsigned int,
//                                           unsigned int> &cell_range) const
//  {
//
//    VectorizedArray<number> sigma_val;
//    unsigned int cell_index;
//    const unsigned int vectorization_length =
//        VectorizedArray<number>::n_array_elements;
//
//    // Only first two template parameters always required
//    // <dimension, fe_degree, quad_degree, n_components, number_type>
//    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> fe_eval (data);
//
//    // Loop over the cells range we are working in this kernel
//    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
//    {
//      cell_index = cell * vectorization_length;
//      // Re-init Values
//      fe_eval.reinit (cell);
//
//      for (unsigned int v = 0; v < data.n_components_filled (cell); ++v)
//        sigma_val[v] = get_coeff(cell_index + v);
//
//      // Read in the values of the source vector including the resolution constraints.
//      // This stores u_cell.
//      fe_eval.read_dof_values (src);
//
//      // Compute unit cell values and gradients (not hessians)
//      fe_eval.evaluate (true, false, false);
//
//      // Jacobi transformation:
//      // Sigma * gradients(in real space) * JxW
//      for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
//      {
//        fe_eval.submit_value (sigma_val * fe_eval.get_value (q), q);
//      }
//      // Summation over all quadrature points in the cell.
//      fe_eval.integrate (true, false);
//
//      // Distribute the local values in the global vector
//      fe_eval.distribute_local_to_global (dst);
//    }
//  }
//
//  template <int dim>
//  template <int fe_degree, typename number>
//  void
//  ScatteringMatrixFree<dim>::vmult_add (Vector<double> &dst,
//                                        const Vector<double> &src) const
//  {
//    data.cell_loop (&ScatteringMatrixFree::local_apply, this, dst, src);
//  }

  template class ScatteringMatrixFree<1> ;
  template class ScatteringMatrixFree<2> ;
  template class ScatteringMatrixFree<3> ;

} // end of namespace Forest

