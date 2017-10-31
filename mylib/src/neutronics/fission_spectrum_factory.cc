/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_factory.cc
 * @brief  Implementation of class template FissionSpectrumFactory
 */

#include "neutronics/fission_spectrum_factory.h"

#include "neutronics/fission_spectrum_matrix_built.h"
#include "neutronics/fission_spectrum_matrix_free.h"

namespace Forest
{
  template <int dim>
  FissionSpectrumPtr<dim>
  FissionSpectrumFactory<dim>::New (State<dim> &state,
                                    const bool matrix_free)
  {
    if (matrix_free)
    {
      return FissionSpectrumPtr<dim> (new FissionSpectrumMatrixFree<dim> (state));
    }
    else
    {
      return FissionSpectrumPtr<dim> (new FissionSpectrumMatrixBuilt<dim> (state));
    }
  }

  template class FissionSpectrumFactory<1> ;
  template class FissionSpectrumFactory<2> ;
  template class FissionSpectrumFactory<3> ;

} // end of namespace Forest
