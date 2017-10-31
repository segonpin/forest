/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op.cc
 * @brief  Implementation of class template TransOp
 */

#include "neutronics/trans_op.h"

#include "algebra/sn_vector.h"        // for SnVector

#include "deal.II/lac/block_vector.h" // for BlockVector

namespace Forest
{
  template <int dim>
  void
  TransOp<dim>::apply_L0_inv (SnVector & phi, const SnVector & source)
  {
    for (unsigned int g = 0; g < get_n_groups (); ++g)
    {
      apply_L0_inv_g (phi.m_data[g], source.m_data[g], g);
    }
  }

  template <int dim>
  void
  TransOp<dim>::apply_L0_inv (SnVector & phi,
                              const BlockVector<double> & source,
                              const unsigned int g)
  {
    apply_L0_inv_g (phi.m_data[g], source, g);
  }

  template <int dim>
  void
  TransOp<dim>::apply_L0_inv (SnVector & phi,
                              const Vector<double> & source,
                              const unsigned int g,
                              const unsigned int i_ord)
  {
    apply_L0_inv_g_iord (phi.m_data[g].block (i_ord), source, g, i_ord);
  }

  template <int dim>
  void
  TransOp<dim>::apply_L0_inv_g (BlockVector<double> & phi,
                                const BlockVector<double> & source,
                                const unsigned int g)
  {
    for (unsigned int i_ord = 0; i_ord < get_n_angles (); ++i_ord)
    {
      apply_L0_inv_g_iord (phi.block (i_ord), source.block (i_ord), g, i_ord);
    }
  }

  template class TransOp<1> ;
  template class TransOp<2> ;
  template class TransOp<3> ;

} // end of namespace Forest
