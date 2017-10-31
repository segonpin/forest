/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_factory.cc
 * @brief  Implementation of class template DiffFactory
 */

#include "neutronics/diff_factory.h"

#include "neutronics/diff_op_matrix_free.h"
#include "neutronics/diff_op_matrix_built.h"

namespace Forest
{
  template <int dim>
  DiffOpPtr<dim>
  DiffFactory<dim>::New (State<dim> &state,
                         const bool matrix_free)
  {
    if (matrix_free)
    {
      return DiffOpPtr<dim> (new DiffOpMatrixFree<dim> (state));
    }
    else
    {
      return DiffOpPtr<dim> (new DiffOpMatrixBuilt<dim> (state));
    }
  }

  template class DiffFactory<1> ;
  template class DiffFactory<2> ;
  template class DiffFactory<3> ;

} // end of namespace Forest
