/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eigen_standard.cc
 * @brief  Implementation of class template EigenStandard
 */

#include "algebra/eigen_standard.h"

#include "algebra/equation.h"
#include "neutronics/state.h"

#include <iomanip>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  void
  EigenStandard<dim>::vmult (SnVector & dst,
                             const SnVector & src) const
  {
    mp_equation->calculate_fission_density (m_fission_density_tmp, src);
    mp_equation->distribute_fission_rate (m_fission_source_tmp,
        m_fission_density_tmp);
    m_solver.iterate (mp_equation, m_phi_tmp, m_fission_source_tmp);
    dst = m_phi_tmp;
  }

  template <int dim>
  void
  EigenStandard<dim>::set_solution_sizes (SnVector & vec) const
  {
    const unsigned int n_groups = mp_state.get_n_groups ();
    vec.reinit (n_groups, mp_equation->get_flux_bv_sizes ());
  }

  template <int dim>
  void
  EigenStandard<dim>::set_initial_guess () const
  {
    // we set the initial value for the flux constantly equal to one
    m_phi_tmp = 1.0;
    /*// we now calculate the norm of the initial guess
     const double scaling_factor = 1./m_phi_tmp.l2_norm();
     // and normalize the previous vectors
     m_phi_tmp *= scaling_factor;*/
  }

  template <int dim>
  SnVector &
  EigenStandard<dim>::get_initial_guess () const
  {
    return m_phi_tmp;
  }

  template <int dim>
  void
  EigenStandard<dim>::update_vectors(const SnVector & src,
      const double ) const
  {
    m_phi_tmp = src;
    mp_equation->calculate_fission_density (m_fission_density_tmp, m_phi_tmp);
    mp_equation->distribute_fission_rate (m_fission_source_tmp,
        m_fission_density_tmp);
  }

  template <int dim>
  double
  EigenStandard<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption +=     m_phi_tmp.memory_consumption() +
        m_fission_density_tmp.memory_consumption() +
        m_fission_source_tmp.memory_consumption() +
        m_solver.memory_consumption();
    return memory_consumption;
  }

  template class EigenStandard<1> ;
  template class EigenStandard<2> ;
  template class EigenStandard<3> ;

} // namespace Forest
