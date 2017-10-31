/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_factory.h
 * @brief  FissionSourceFactory class template declarations
 */

#ifndef FOREST_FISSION_SOURCE_FACTORY_H
#define FOREST_FISSION_SOURCE_FACTORY_H

#include "neutronics/fission_source.h"
#include "neutronics/state_fwd.h"

#include <memory>

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using FissionSourcePtr = std::shared_ptr<FissionSource<dim> >;

  /**
   * @class FissionSourceFactory
   * @brief Factory to create fission source operators (w/o matrix)
   */
  template <int dim>
  class FissionSourceFactory
  {
  public:
    /** @todo document me */
    static FissionSourcePtr<dim>
    New (State<dim> &state,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SOURCE_FACTORY_H */
