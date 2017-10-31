/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eigen_prob_factory.h
 * @brief   EigenProbFactory class template declarations
 */

#ifndef FOREST_EIGEN_PROB_FACTORY_H
#define FOREST_EIGEN_PROB_FACTORY_H

#include "neutronics/state_fwd.h"       // For class State&
#include "algebra/equation_fwd.h"       // for Equation &

#include <memory>

#ifndef DOXYGEN
namespace Forest { template <int dim> class EigenProb; }
#endif

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using EigenProbPtr = std::shared_ptr<EigenProb<dim> >;

  /**
   * @ingroup ForestAlgebra
   * @class EigenProbFactory
   * @brief Factory class creating EigenProbPtr
   */
  template <int dim>
  class EigenProbFactory
  {
  public:

    /**
     * @brief We build the Eigenvalue problem
     * @details We choose between Standard or Fission Density problem
     * for the moment.
     */
    static EigenProbPtr<dim>
    New (State<dim> &state,
         std::shared_ptr<Equation<dim> > & equation,
         const bool use_fission_density = true);
  };
} /* end of namespace Forest */

#endif /* FOREST_EIGEN_PROB_FACTORY_H */
