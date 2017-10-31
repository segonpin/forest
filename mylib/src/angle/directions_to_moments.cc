/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/directions_to_moments.cc
 * @brief  Implementation of class template DirectionsToMoments
 */

#include "angle/directions_to_moments.h"

#include "angle/quadrature_base.h"
#include "angle/spherical_harmonics.h"

#include "deal.II/lac/vector.h"         // for Vector
#include "deal.II/lac/block_vector.h"   // for BlockVector

#include <cmath>                        // for sqrt
#include <vector>                       // for vector
#include <memory>                       // for shared_ptr

namespace Forest
{
  using namespace dealii;

  template <int dim>
  DirectionsToMoments<dim>::DirectionsToMoments (std::shared_ptr<
                                                 QuadratureBase<dim> > & ord,
                                                 unsigned int moments_order)
      : mp_ord (ord),
        m_n_angles (mp_ord->get_n_angles ()),
        m_moments_order (moments_order),
        m_moment_indexer (m_moments_order),
        m_n_moments (m_moment_indexer.get_number_of_moments ())
  {
  }

  template <int dim>
  void
  DirectionsToMoments<dim>::vmult (BlockVector<double> & phi,
                                   const BlockVector<double> & psi)
  {
    for (unsigned int i = 0; i < m_n_moments; ++i)
    {
      vmult (phi.block (i), psi, i);
    }
  }

  template <int dim>
  void
  DirectionsToMoments<dim>::vmult (Vector<double> & phi_i,
                                   const BlockVector<double> & psi,
                                   const unsigned int i)
  {
    phi_i = 0;
    for (unsigned int n = 0; n < m_n_angles; ++n)
    {
      vmult_add (phi_i, psi.block (n), i, n);
    }
  }

  template <int dim>
  void
  DirectionsToMoments<dim>::vmult_add (BlockVector<double> & phi,
                                       const Vector<double> & psi_n,
                                       const unsigned int n)
  {
    for (unsigned int i = 0; i < m_n_moments; ++i)
    {
      vmult_add (phi.block (i), psi_n, i, n);
    }
  }

  template <int dim>
  void
  DirectionsToMoments<dim>::vmult_add (Vector<double> & phi_i,
                                       const Vector<double> & psi_n,
                                       const unsigned int i,
                                       const unsigned int n)
  {
    std::vector<int> lm = {m_moment_indexer.get_l (i), m_moment_indexer.get_m (i)};
    std::vector<double> q = mp_ord->get_q (n);
    q[2] = std::sqrt (1 - (q[0] * q[0] + q[1] * q[1]));
    const double w_n = mp_ord->get_w (n);
    const double y_in = SphericalHarmonics<dim>::Y_lm (lm, q);
    phi_i.add (y_in*w_n, psi_n);
  }

  template class DirectionsToMoments<1> ;
  template class DirectionsToMoments<2> ;
  template class DirectionsToMoments<3> ;

} // end of namespace Forest
