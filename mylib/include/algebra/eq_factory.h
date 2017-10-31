/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eq_factory.h
 * @brief   EqFactory class template declarations
 */

#ifndef FOREST_EQ_FACTORY_H
#define FOREST_EQ_FACTORY_H

#include "algebra/equation_fwd.h"       // for Equation &
#include "neutronics/state_fwd.h"       // For class State&

#include <memory>

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using EqPtr = std::shared_ptr<Equation<dim> >;

  /**
   * @brief Factory to create equations (transport/diffusion, with/without matrix)
   * @ingroup ForestAlgebra
   */
  template <int dim>
  class EqFactory
  {
  public:

    /** @todo document me */
    static EqPtr<dim>
    New (State<dim> &state,
         bool transport = true);
  };
} // end of namespace Forest

#endif /* FOREST_EQ_FACTORY_H */
