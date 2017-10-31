/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_factory.h
 * @brief  TransFactory class template declarations
 */

#ifndef FOREST_TRANS_FACTORY_H
#define FOREST_TRANS_FACTORY_H

#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#include <memory>

#ifndef DOXYGEN
namespace Forest { template <int dim> class TransOp; }
#endif

namespace Forest
{

  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using TransOpPtr = std::shared_ptr<TransOp<dim> >;

  /**
   * @class TransFactory
   * @ingroup ForestNeutronics
   * @brief Factory to create Transport operators (with/without matrix)
   */
  template <int dim>
  class TransFactory
  {
  public:

    /** @todo document me */
    static TransOpPtr<dim>
    New (State<dim> &state,
         std::shared_ptr<QuadratureBase<dim> > & ord,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_TRANS_FACTORY_H */
