/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_matrix_built.cc
 * @brief  Implementation of class template FissionSpectrumMatrixBuilt
 */

#include "neutronics/fission_spectrum_matrix_built.h"

#include "neutronics/state.h"
#include "input/input.h"
#include "algebra/sn_vector.h"             // for SnVector

#include "deal.II/dofs/dof_tools.h"      // for FiniteElement, etc
#include "deal.II/dofs/dof_handler.h"              // for DoFHandler

namespace Forest
{
  using namespace dealii;

  template <int dim>
  FissionSpectrumMatrixBuilt<dim>::FissionSpectrumMatrixBuilt (State<dim> & state)
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
  FissionSpectrumMatrixBuilt<dim>::set_non_zero_pattern ()
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

  template <int dim>
  void
  FissionSpectrumMatrixBuilt<dim>::build ()
  {
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_dofs = mp_dof_handler.n_dofs ();
    m_chi_matrix.resize (n_groups, Vector<double> (n_dofs));

    const FiniteElement<dim> & fe (mp_state.get_fe ());
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    double chi;
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
            m_chi_matrix[g][local_dof_indices[i]] = chi;
          }
        }
        else
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            m_chi_matrix[g][local_dof_indices[i]] = 0.0;
          }
        }
      }
    }
  }

  template <int dim>
  void
  FissionSpectrumMatrixBuilt<dim>::vmult (SnVector & source,
                                          const SnVector & phi)
  {
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_dofs = mp_dof_handler.n_dofs ();
    // and here we apply the fission spectrum
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      if (m_non_zero[g])
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
        {
          source.m_data[g].block (0)[i] = m_chi_matrix[g][i]
              * phi.m_data[0].block (0)[i];
        }
      }
      else
      {
        source.m_data[g].block (0) = 0.0;
      }
    }
  }

  template <int dim>
  void
  FissionSpectrumMatrixBuilt<dim>::vmult_add (SnVector & source,
                                              const SnVector & phi)
  {
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_dofs = mp_dof_handler.n_dofs ();
    // and here we apply the fission spectrum
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      if (m_non_zero[g])
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
        {
          source.m_data[g].block (0)[i] += m_chi_matrix[g][i]
              * phi.m_data[0].block (0)[i];
        }
      }
    }
  }

  template <int dim>
  double
  FissionSpectrumMatrixBuilt<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      memory_consumption += m_chi_matrix[g].memory_consumption ();
    }
    return memory_consumption;
  }

  template class FissionSpectrumMatrixBuilt<1> ;
  template class FissionSpectrumMatrixBuilt<2> ;
  template class FissionSpectrumMatrixBuilt<3> ;

} // end of namespace Forest
