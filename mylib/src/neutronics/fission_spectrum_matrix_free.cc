/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_matrix_free.cc
 * @brief  Implementation of class template FissionSpectrumMatrixFree
 */

#include "neutronics/fission_spectrum_matrix_free.h"

#include "neutronics/state.h"
#include "input/input.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/dofs/dof_tools.h"      // for FiniteElement, etc
#include "deal.II/dofs/dof_handler.h"              // for DoFHandler

#include <vector>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  FissionSpectrumMatrixFree<dim>::FissionSpectrumMatrixFree (State<dim> & state)
      : mp_state (state),
        mp_dof_handler (state.get_dof_handler ())
  {
    // Set the fission spectrum non zero pattern
    this->set_non_zero_pattern ();
  }

  template <int dim>
  void
  FissionSpectrumMatrixFree<dim>::set_non_zero_pattern ()
  {
    // we consider zero anything below this value;
    const double eps = 1.e-14;
    const unsigned int n_groups = mp_state.get_n_groups ();
    m_non_zero.resize (n_groups, false);
    double chi;
    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      for (cell = mp_dof_handler.begin_active (); cell != endc; ++cell)
      {
        /*
         * Here the chi should come from TH
         */
        chi = mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_chi[g];
        if (std::abs (chi) > eps)
        {
          m_non_zero[g] = true;
          // Once we find a non-zero value, then stop the loop for this group
          continue;
        }
      }
    }
  }

  /**
   * @todo prepare a function that apply the fission spectrum depending on
   * the material, running over the whole triangulation
   */
  template <int dim>
  void
  FissionSpectrumMatrixFree<dim>::vmult (SnVector & source,
                                         const SnVector & phi)
  {
    // We need these things
    const unsigned int n_groups = mp_state.get_n_groups ();
    const FiniteElement<dim> & fe (mp_state.get_fe ());
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    double chi;
    // running over the cells to multiply by chi
    // container for the cross section and indices
    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (; cell != endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        if (m_non_zero[g])
        {
          /*
           * Here the chi should come from TH
           */
          chi = mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_chi[g];
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int j = local_dof_indices[i];
            source.m_data[g].block (0)[j] = chi * phi.m_data[0].block (0)[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            source.m_data[g].block (0)[local_dof_indices[i]] = 0.0;
          }
        }
      }
    }
  }

  /**
   * @todo prepare a function that apply the fission spectrum depending on
   * the material, running over the whole triangulation
   */
  template <int dim>
  void
  FissionSpectrumMatrixFree<dim>::vmult_add (SnVector & source,
                                             const SnVector & phi)
  {
    // We need these things
    const unsigned int n_groups = mp_state.get_n_groups ();
    const FiniteElement<dim> & fe (mp_state.get_fe ());
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    // container for the cross section and indices
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    double chi;
    // running over the cells to multiply by chi
    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (; cell != endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        if (m_non_zero[g])
        {
          /*
           * Here the chi should come from TH
           */
          chi = mp_state.mp_data.mp_mat.m_xs.at (cell->material_id ()).m_chi[g];
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int j = local_dof_indices[i];
            source.m_data[g].block (0)[j] += chi * phi.m_data[0].block (0)[j];
          }
        }
      }
    }
  }

  template <int dim>
  double
  FissionSpectrumMatrixFree<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    return memory_consumption;
  }

  template class FissionSpectrumMatrixFree<1> ;
  template class FissionSpectrumMatrixFree<2> ;
  template class FissionSpectrumMatrixFree<3> ;

} // end of namespace Forest
