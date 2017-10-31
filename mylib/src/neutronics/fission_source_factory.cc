/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_factory.cc
 * @brief  Implementation of class template FissionSourceFactory
 */

#include "neutronics/fission_source_factory.h"

#include "neutronics/fission_source_matrix_built.h"
#include "neutronics/fission_source_matrix_free.h"

namespace Forest
{
  template <int dim>
  FissionSourcePtr<dim>
  FissionSourceFactory<dim>::New (State<dim> &state,
                                  const bool matrix_free)
  {
    if (matrix_free)
    {
      return FissionSourcePtr<dim> (new FissionSourceMatrixFree<dim> (state));
    }
    else
    {
      return FissionSourcePtr<dim> (new FissionSourceMatrixBuilt<dim> (state));
    }
  }

  template class FissionSourceFactory<1> ;
  template class FissionSourceFactory<2> ;
  template class FissionSourceFactory<3> ;

} // end of namespace Forest
