/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/moment_indexer.cc
 * @brief  Implementation of class template MomentIndexer
 */

#include "angle/moment_indexer.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <iostream>
#include <vector>
#include <memory>

namespace Forest
{

  template <int dim>
  MomentIndexer<dim>::MomentIndexer (const unsigned int angular_order)
      : m_angular_order (angular_order),
        m_number_of_moments(0)
  {
    this->construct ();

    // Fill the l and m vectors
    m_l_vec.resize (m_number_of_moments, 0);
    m_m_vec.resize (m_number_of_moments, 0);
    unsigned int j = 0;
    for (unsigned int l = 0; l <= m_angular_order; ++l)
    {
      for (unsigned int i = 0; i < m_m_index[l].size (); ++i, ++j)
      {
        m_l_vec[j] = l;
        m_m_vec[j] = m_m_index[l][i];
      }
    }
  }

  template <int dim>
  void
  MomentIndexer<dim>::display () const
  {
    std::cout << " MOMENTS INDEXER: " << std::endl;
    for (unsigned int i = 0; i < m_number_of_moments; ++i)
    {
      const int l = this->get_l(i);
      const int m = this->get_m(i);
      std::cout << " " << i << " " << l << " " << m << std::endl;
    }
  }

  template <int dim>
  void
  MomentIndexer<dim>::construct ()
  {
    Assert(false, dealii::ExcMessage("Wrong dimension provided."));
  }

  template <>
  void
  MomentIndexer<1>::construct ()
  {
    m_number_of_moments = m_angular_order + 1;
    m_m_index.resize (m_angular_order + 1);
    for (unsigned int l = 0; l <= m_angular_order; ++l)
    {
      m_m_index[l].resize (1, 0); /* m is always zero in 1D */
    }
  }

  template <>
  void
  MomentIndexer<2>::construct ()
  {
    m_number_of_moments = (m_angular_order + 1) * (m_angular_order + 2) / 2;
    m_m_index.resize (m_angular_order + 1);
    for (int l = 0; l <= (int) m_angular_order; ++l)
    {
      for (int m = -l; m <= l; m++)
      {
        if ((m + l) % 2 == 0)
        {
          m_m_index[l].push_back (m);
        }
      }
    }
  }

  template <>
  void
  MomentIndexer<3>::construct ()
  {
    m_number_of_moments = (m_angular_order + 1) * (m_angular_order + 1);
    m_m_index.resize (m_angular_order + 1);
    for (int l = 0; l <= (int) m_angular_order; ++l)
    {
      for (int m = -l; m <= l; ++m)
      {
        m_m_index[l].push_back (m);
      }
    }
  }

  /* Explicit instantiations. */
  template class MomentIndexer<1> ;
  template class MomentIndexer<2> ;
  template class MomentIndexer<3> ;
}
// end of namespace Forest
