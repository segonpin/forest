/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_matrix_free.cc
 * @brief  Implementation of class template FissionSourceMatrixFree
 */

#include "neutronics/fission_source_matrix_free.h"

#include "neutronics/state.h"
#include "input/input.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "deal.II/fe/fe.h"               // for FiniteElement
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc
#include "deal.II/fe/fe_values.h"        // for FEValues
#include "deal.II/grid/tria.h"           // for Triangulation
#include "deal.II/grid/grid_generator.h" // for GridGenerator
#include "deal.II/lac/vector.h"

#include <vector>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  FissionSourceMatrixFree<dim>::FissionSourceMatrixFree (State<dim> & state)
      : mp_state (state),
        mp_dof_handler (state.get_dof_handler ())
  {
    // We first setup the non zero pattern, in order to avoid operations later
    this->set_non_zero_pattern ();
  }

  template <int dim>
  void
  FissionSourceMatrixFree<dim>::set_non_zero_pattern ()
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
  FissionSourceMatrixFree<dim>::vmult (SnVector & source,
                                       const SnVector & phi)
  {
    source.m_data[0].block (0) = 0;
    vmult_add (source, phi);
  }

  template <int dim>
  void
  FissionSourceMatrixFree<dim>::vmult_add (SnVector & source,
                                           const SnVector & phi)
  {
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      if (m_non_zero[g])
      {
        ingroup_fission_add (source.m_data[0].block (0),
            phi.m_data[g].block (0), g);
      }
    }
  }

  template <int dim>
  void
  FissionSourceMatrixFree<dim>::ingroup_fission_add (Vector<double> & source,
                                                     const Vector<double> & phi,
                                                     const unsigned int g)
  {
    const FiniteElement<dim> & fe (mp_state.get_fe ());
    QGauss<dim> quadrature_formula (fe.degree + 1);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size ();
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell);
    Vector<double> temp_vector (n_q_points);

    FEValues<dim> fe_values_ref (fe, quadrature_formula, update_values);
    Triangulation<dim, dim> reference_cell;
    GridGenerator::hyper_cube (reference_cell, 0., 1.);
    fe_values_ref.reinit (reference_cell.begin ());
    FEValues<dim> fe_values (fe, quadrature_formula,
        update_JxW_values | update_quadrature_points);

    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      /*
       * Here the nusigf should come from TH
       */
      double nusigf =
          mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_nusigf[g];

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        cell_src (i) = phi (local_dof_indices[i]);
      }
      temp_vector = 0;
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          temp_vector (q) += fe_values_ref.shape_value (i, q) * cell_src (i);
        }
      }

      // multiply by coefficient and integration weight
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        temp_vector (q) *= fe_values.JxW (q) * nusigf;
      }

      cell_dst = 0;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          cell_dst (i) += fe_values_ref.shape_value (i, q) * temp_vector (q);
        }
      }
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        source (local_dof_indices[i]) += cell_dst (i);
      }
    }
  }

  template <int dim>
  double
  FissionSourceMatrixFree<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    return memory_consumption;
  }

  template class FissionSourceMatrixFree<1> ;
  template class FissionSourceMatrixFree<2> ;
  template class FissionSourceMatrixFree<3> ;

} // end of namespace Forest
