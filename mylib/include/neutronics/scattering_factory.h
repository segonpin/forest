/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_factory.h
 * @brief  ScatteringFactory class template declarations
 */

#ifndef FOREST_SCATTERING_FACTORY_H
#define FOREST_SCATTERING_FACTORY_H

#include "neutronics/scattering_source.h"

#include "neutronics/state_fwd.h"       // For class State&

#include <memory>

namespace Forest
{

  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using ScatteringSourcePtr = std::shared_ptr<ScatteringSource<dim> >;

  /**
   * @class ScatteringFactory
   * @ingroup ForestNeutronics
   * @brief Factory to create Scattering Source operators (with/without matrix)
   */
  template <int dim>
  class ScatteringFactory
  {
  public:
    /** @todo document me */
    static ScatteringSourcePtr<dim>
    New (State<dim> &state,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_SCATTERING_FACTORY_H */
