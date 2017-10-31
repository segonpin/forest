/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/boundary_values_factory.h
 * @brief  BoundaryValuesFactory class template declarations
 */

#ifndef FOREST_BOUNDARY_VALUES_FACTORY_H
#define FOREST_BOUNDARY_VALUES_FACTORY_H

#include "neutronics/boundaryvalues.h"
#include "neutronics/state_fwd.h"
#include "angle/quadrature_base.h"

#include <memory>

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using BoundaryValuesPtr = std::shared_ptr<BoundaryValues<dim> >;

  /**
   * @ingroup NeutronTransport
   * @class BoundaryValuesFactory
   * @brief Factory to create boundary values (w/o matrix)
   */
  template <int dim>
  class BoundaryValuesFactory
  {
  public:
    /** @todo document me */
    static BoundaryValuesPtr<dim>
    New (State<dim> & state,
         std::shared_ptr<QuadratureBase<dim> > & ord,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_BOUNDARY_VALUES_FACTORY_H */
