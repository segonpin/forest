/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_factory.cc
 * @brief  Implementation of class template TransFactory
 */

#include "neutronics/trans_factory.h"

#include "neutronics/trans_op_matrix_built.h"
#include "neutronics/trans_op_matrix_free.h"
#include "neutronics/trans_op_mf_pre.h"

#include "geometry/prob_geom.h"
#include "neutronics/state.h"

namespace Forest
{
  template <int dim>
  TransOpPtr<dim>
  TransFactory<dim>::New (State<dim> &state,
      std::shared_ptr<QuadratureBase<dim> > &ord,
      const bool matrix_free)
  {
    if (matrix_free)
    {
      if (state.mp_geom.get_mesh_regularity() == 0)
      {
        return TransOpPtr<dim>(new TransOpMatrixFree<dim> (state,ord));
      }
      else
      {
        return TransOpPtr<dim>(new TransOpMatrixFreePre<dim> (state,ord));
      }
    }
    else
    {
      return TransOpPtr<dim>(new TransOpMatrixBuilt<dim> (state,ord));
    }
  }

  template class TransFactory<1> ;
  template class TransFactory<2> ;
  template class TransFactory<3> ;

} // end of namespace Forest
