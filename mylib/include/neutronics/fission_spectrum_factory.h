/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_factory.h
 * @brief  FissionSpectrumFactory class template declarations
 */

#ifndef FOREST_FISSION_SPECTRUM_FACTORY_H
#define FOREST_FISSION_SPECTRUM_FACTORY_H

#include "neutronics/fission_spectrum.h"
#include "neutronics/state_fwd.h"

#include <memory>

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using FissionSpectrumPtr = std::shared_ptr<FissionSpectrum<dim> >;

  /**
   * @class FissionSpectrumFactory
   * @brief Factory to create fission spectrum operators (w/o matrix)
   */
  template <int dim>
  class FissionSpectrumFactory
  {
  public:
    /** @todo document me */
    static FissionSpectrumPtr<dim>
    New (State<dim> &state,
         bool matrix_free = true);
  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SPECTRUM_FACTORY_H */
