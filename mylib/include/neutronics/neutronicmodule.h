/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/neutronicmodule.h
 * @brief  NeutronicModule class template declarations
 */

#ifndef FOREST_NEUTRONICMODULE_H
#define FOREST_NEUTRONICMODULE_H

#include "neutronics/state_fwd.h"       // For class State&
//#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&
//#include "angle/quadrature_base_fwd.h"

#ifndef DOXYGEN
namespace Forest { template <int dim> class EigenSolver; }
#endif

namespace Forest
{

  //---------------------------------------------------------------------------
  // Class using SLEPc to solve the problem
  //---------------------------------------------------------------------------

  /**
   * @class NeutronicModule
   * @ingroup ForestNeutronics
   * @brief A Class related to neutron transport operators and solvers
   *
   * Here we solve the Neutron Transport equation
   *
   * \f{align*}{
   * \frac{1}{v(E)} \frac{\partial}{\partial t} \Psi(\mathbf{r},\Omega,E,t) +
   * \Omega \cdot \nabla \Psi(\mathbf{r},\Omega,E,t) +
   * \Sigma_{T}(\mathbf{r},E,t) \Psi(\mathbf{r},\Omega,E,t) = \nonumber \\
   * \int_{(4\pi)} \int_0^{\infty}
   * \Sigma_{s}(\mathbf{r},\Omega^{\prime}\to\Omega,E^{\prime}\to E,t)
   * \Psi(\mathbf{r},\Omega^{\prime},E^{\prime},t)
   * d\Omega^{\prime} dE^{\prime} + \nonumber \\
   * \frac{\chi(E)}{4\pi}\int_{0}^{\infty}
   * \nu(E^\prime)\Sigma_{f}(\mathbf{r},E^{\prime},t)
   * \Phi(\mathbf{r},E^{\prime},t) dE^{\prime} ,
   * \f}
   *
   * or any approximation of it. Different methods can be used for the
   * angular discretization, like for example the discrete ordinates,
   * \f$ S_N \f$, or the Spherical Harmonics method, \f$ P_N \f$. Similarly,
   * different methods can be used for the spatial approximation, as the
   * Finite Element Method here.
   *
   * In order to solve the equation, we need the following
   *  - State object
   *  - Problem Data
   *  - Geometry
   *  - Quadrature (angular discretization)
   *  - Method object (right now it is the neutron transport solver)
   *
   * Once we have all the previous components, we are ready to solve the
   * problem.
   *
   * @todo The Quadrature should be created inside, should not be coming
   * from outside. The reason to not take it from outside is because it
   * will not be always needed, but when needed it is purely neutronic.
   *
   */

  template <int dim>
  class NeutronicModule
  {
  public:

    /** @todo document me */
    NeutronicModule (State<dim> & state)
        : mp_state (state)
    {
    }

    /** @todo document me */
    void
    run ();

  private:

    State<dim> & mp_state;
    // QuadratureBase<dim> m_ord;
    // EigenSolver<dim> m_problem;
  };

} // end of namespace Forest
#endif
