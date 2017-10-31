/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eq_factory.cc
 * @brief  Implementation of class template EqFactory
 */

#include "algebra/eq_factory.h"
#include "algebra/eq_transport.h"
#include "algebra/eq_diffusion.h"

namespace Forest
{
  template <int dim>
  EqPtr<dim>
  EqFactory<dim>::New (State<dim> &state,
                       const bool transport)
  {
    if (transport)
      return EqPtr<dim> (new EqTransport<dim> (state));
    else
      return EqPtr<dim> (new EqDiffusion<dim> (state));
  }
  template class EqFactory<1> ;
  template class EqFactory<2> ;
  template class EqFactory<3> ;
} // end of namespace Forest
