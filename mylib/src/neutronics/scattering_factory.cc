/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_factory.cc
 * @brief  Implementation of class template ScatteringFactory
 */

#include "neutronics/scattering_factory.h"

#include "neutronics/scattering_matrix_built.h"
#include "neutronics/scattering_matrix_free.h"

namespace Forest
{
  template <int dim>
  ScatteringSourcePtr<dim>
  ScatteringFactory<dim>::New (State<dim> &state,
                               const bool matrix_free)
  {
    if (matrix_free)
    {
      return ScatteringSourcePtr<dim> (new ScatteringMatrixFree<dim> (state));
    }
    else
    {
      return ScatteringSourcePtr<dim> (new ScatteringMatrixBuilt<dim> (state));
    }
  }

  template class ScatteringFactory<1> ;
  template class ScatteringFactory<2> ;
  template class ScatteringFactory<3> ;

} // end of namespace Forest
