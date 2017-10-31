/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eigen_fission_density.cc
 * @brief  Implementation of class template EigenFissionDensity
 */

#include "algebra/equation.h"
#include "algebra/eigen_fission_density.h"
#include "neutronics/state.h"

#include <iomanip>
#include <memory>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  void
  EigenFissionDensity<dim>::vmult (SnVector & dst,
                                   const SnVector & src) const
  {
    /* We distribute the fission rate using the fission spectra. */
    mp_equation->distribute_fission_rate (m_fission_source_tmp, src);
    /* Solve the source problem with the source iteration. */
    m_solver.iterate (mp_equation, m_phi_tmp, m_fission_source_tmp);
    /* We generate the new fission rate. */
    mp_equation->calculate_fission_density (m_fission_density_tmp, m_phi_tmp);
    dst = m_fission_density_tmp;
  }

  template <int dim>
  void
  EigenFissionDensity<dim>::set_solution_sizes (SnVector & vec) const
  {
    vec.reinit (1, mp_equation->get_fission_density_sizes ());
  }

  template <int dim>
  void
  EigenFissionDensity<dim>::set_initial_guess () const
  {
    // we set the initial value for the flux constantly equal to one
    m_phi_tmp = 1.0;
    // and then we calculate the fission density
    mp_equation->calculate_fission_density (m_fission_density_tmp, m_phi_tmp);
    // we now calculate the norm of the initial guess
    const double scaling_factor = 1. / m_fission_density_tmp.l2_norm ();
    // and normalize the previous vectors
    m_fission_density_tmp *= scaling_factor;
    m_phi_tmp *= scaling_factor;
  }

  template <int dim>
  SnVector &
  EigenFissionDensity<dim>::get_initial_guess () const
  {
    return m_fission_density_tmp;
  }

  template <int dim>
  void
  EigenFissionDensity<dim>::update_vectors(const SnVector & src,
      const double keff) const
  {
    m_fission_density_tmp = src;
    mp_equation->distribute_fission_rate (m_fission_source_tmp, src);
    m_solver.iterate (mp_equation, m_phi_tmp, m_fission_source_tmp);
    m_phi_tmp *= 1.0/keff;
  }

  template <int dim>
  double
  EigenFissionDensity<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption +=     m_phi_tmp.memory_consumption() +
        m_fission_density_tmp.memory_consumption() +
        m_fission_source_tmp.memory_consumption() +
        m_solver.memory_consumption();
    return memory_consumption;
  }

  template class EigenFissionDensity<1> ;
  template class EigenFissionDensity<2> ;
  template class EigenFissionDensity<3> ;

} /* namespace Forest */
