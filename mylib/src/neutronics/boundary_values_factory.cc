/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/boundary_values_factory.cc
 * @brief  Implementation of class template BoundaryValuesFactory
 */

#include "neutronics/boundary_values_factory.h"

#include "neutronics/boundaryvalues.h"

namespace Forest
{
  template <int dim>
  BoundaryValuesPtr<dim>
  BoundaryValuesFactory<dim>::New (State<dim> & state,
                                   std::shared_ptr<QuadratureBase<dim> > & ord,
                                   const bool matrix_free)
  {
    if (matrix_free)
    {
      return BoundaryValuesPtr<dim> (new BoundaryValues<dim> (state, ord));
    }
    else
    {
      return BoundaryValuesPtr<dim> (new BoundaryValues<dim> (state, ord));
    }
  }

  template class BoundaryValuesFactory<1> ;
  template class BoundaryValuesFactory<2> ;
  template class BoundaryValuesFactory<3> ;

} // end of namespace Forest
