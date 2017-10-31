/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/moments_to_directions.cc
 * @brief  Implementation of class template MomentsToDirections
 */

#include "angle/moments_to_directions.h"

#include "angle/quadrature_base.h"
#include "angle/spherical_harmonics.h"

#include "deal.II/lac/vector.h"         // for Vector
#include "deal.II/lac/block_vector.h"   // for BlockVector

#include <vector>
#include <cmath>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  MomentsToDirections<dim>::MomentsToDirections (const std::shared_ptr<
                                                     QuadratureBase<dim> > & ord,
                                                 const unsigned int moments_order)
      : mp_ord (ord),
        m_n_angles (mp_ord->get_n_angles ()),
        m_moments_order (moments_order),
        m_moment_indexer (m_moments_order),
        m_n_moments (m_moment_indexer.get_number_of_moments ())
  {
  }

  template <int dim>
  void
  MomentsToDirections<dim>::vmult (BlockVector<double> & psi,
                                   const BlockVector<double> & phi)
  {
    for (unsigned int n = 0; n < m_n_angles; ++n)
    {
      vmult (psi.block (n), phi, n);
    }
  }

  template <int dim>
  void
  MomentsToDirections<dim>::vmult (Vector<double> & psi_n,
                                   const BlockVector<double> & phi,
                                   const unsigned int n)
  {
    this->vmult0 (psi_n, phi.block (0));
    for (unsigned int i = 1; i < m_n_moments; ++i)
    {
      vmult_add (psi_n, phi.block (i), n, i);
    }
  }

  template <int dim>
  void
  MomentsToDirections<dim>::vmult_add (Vector<double> & psi_n,
                                       const Vector<double> & phi_i,
                                       const unsigned int n,
                                       const unsigned int i)
  {
    std::vector<int> lm =
    { m_moment_indexer.get_l (i), m_moment_indexer.get_m (i) };
    std::vector<double> q = mp_ord->get_q (n);
    q[2] = std::sqrt (1 - (q[0] * q[0] + q[1] * q[1]));
    const double y_in = SphericalHarmonics<dim>::Y_lm (lm, q);
    psi_n.add (y_in, phi_i);
  }

  template <int dim>
  void
  MomentsToDirections<dim>::vmult0 (Vector<double> & psi_n,
                                    const Vector<double> & phi_i)
  {
    psi_n = phi_i;
  }

  template class MomentsToDirections<1> ;
  template class MomentsToDirections<2> ;
  template class MomentsToDirections<3> ;

} // end of namespace Forest
