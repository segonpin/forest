/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_factory.h
 * @brief  DiffFactory class template declarations
 */

#ifndef FOREST_DIFF_FACTORY_H
#define FOREST_DIFF_FACTORY_H

#include "neutronics/state_fwd.h"       // For class State&

#include <memory>

#ifndef DOXYGEN
namespace Forest { template <int dim> class DiffOp; }
#endif

namespace Forest
{

  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using DiffOpPtr = std::shared_ptr<DiffOp<dim> >;

  /**
   * @class DiffFactory
   * @ingroup ForestNeutronics
   * @brief Factory to create Diffusion operator (with/without matrix)
   */
  template <int dim>
  class DiffFactory
  {
  public:

    /** @todo document me */
    static DiffOpPtr<dim>
    New (State<dim> &state,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_DIFF_FACTORY_H */
