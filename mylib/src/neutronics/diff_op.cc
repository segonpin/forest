/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op.cc
 * @brief  Implementation of class template DiffOp
 */

#include "neutronics/diff_op.h"

#include "algebra/sn_vector.h"        // for SnVector

#include "deal.II/lac/block_vector.h" // for BlockVector
#include "deal.II/lac/vector.h"       // for Vector

namespace Forest
{
  template <int dim>
  void
  DiffOp<dim>::vmult (SnVector & phi,
                      const SnVector & source)
  {
    for (unsigned int g = 0; g < get_n_groups (); ++g)
    {
      vmult (phi.m_data[g], source.m_data[g], g);
    }
  }

  template <int dim>
  void
  DiffOp<dim>::vmult (SnVector & phi,
                      const BlockVector<double> & source,
                      const unsigned int g)
  {
    vmult (phi.m_data[g], source, g);
  }

  template class DiffOp<1> ;
  template class DiffOp<2> ;
  template class DiffOp<3> ;

} // end of namespace Forest
