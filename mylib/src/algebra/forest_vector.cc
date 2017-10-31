/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/forest_vector.cc
 * @brief  Implementation of classes ForestVector (and MyBlockMatrix)
 */

#include "algebra/forest_vector.h"

#include "deal.II/base/exceptions.h"    // for AssertDimension, Assert

#include <cmath>
#include <cassert>

namespace Forest
{
  using namespace dealii;

  template <typename Number>
  ForestVector<Number>::ForestVector ()
      : m_initialized (false),
        m_n_multi_blocks (0),
        m_block_sizes (0, 0),
        m_data (m_n_multi_blocks, BlockVector<Number> (m_block_sizes))
  {
  }

  template <typename Number>
  ForestVector<Number>::ForestVector (const unsigned int n_multi_blocks,
                                      const std::vector<unsigned int> & block_sizes)
      : m_initialized (true),
        m_n_multi_blocks (n_multi_blocks),
        m_block_sizes (block_sizes),
        m_data (m_n_multi_blocks, BlockVector<Number> (m_block_sizes))
  {
  }

  template <typename Number>
  void
  ForestVector<Number>::reinit (const unsigned int n_multi_blocks,
                                const std::vector<unsigned int> & block_sizes,
                                bool leave_elements_uninitialized)
  {
    m_initialized = true;
    m_n_multi_blocks = n_multi_blocks;
    m_block_sizes = block_sizes;

    m_data.resize (m_n_multi_blocks);
    for (unsigned int i = 0; i < m_n_multi_blocks; ++i)
    {
      m_data[i].reinit (m_block_sizes, leave_elements_uninitialized);
    }
  }

  template <typename Number>
  void
  ForestVector<Number>::reinit (const ForestVector<Number> & v,
                                const bool leave_elements_uninitialized)
  {
    m_initialized = true;
    m_n_multi_blocks = v.m_n_multi_blocks;
    m_block_sizes = v.m_block_sizes;

    m_data.resize (v.m_data.size ());
    for (unsigned int i = 0; i < v.m_data.size (); ++i)
    {
      m_data[i].reinit (v.m_data[i], leave_elements_uninitialized);
    }
  }

  template <typename Number>
  Number
  ForestVector<Number>::operator * (const ForestVector<Number> &v) const
  {
    assert(m_initialized);
    Number val = 0.;
    AssertDimension(size0(), v.size0());
    for (unsigned int i = 0; i < v.size0 (); ++i)
    {
      val += m_data[i] * v.m_data[i];
    }
    return val;
  }

  template <typename Number>
  void
  ForestVector<Number>::add (const ForestVector<Number> &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    for (unsigned int i = 0; i < x.size0 (); ++i)
    {
      m_data[i] += x.m_data[i];
    }
  }

  template <typename Number>
  void
  ForestVector<Number>::add (const Number a,
                             const ForestVector<Number> &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    for (unsigned int i = 0; i < x.size0 (); ++i)
    {
      m_data[i].add (a, x.m_data[i]);
    }
  }

  template <typename Number>
  void
  ForestVector<Number>::sadd (const Number a,
                              const Number b,
                              const ForestVector<Number> &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    for (unsigned int i = 0; i < x.size0 (); ++i)
    {
      m_data[i].sadd (a, b, x.m_data[i]);
    }
  }

  template <typename Number>
  void
  ForestVector<Number>::equ (const Number a,
                             const ForestVector<Number> &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    for (unsigned int i = 0; i < x.size0 (); ++i)
    {
      m_data[i].equ (a, x.m_data[i]);
    }
  }

  template <typename Number>
  ForestVector<Number> &
  ForestVector<Number>::operator *= (const Number a)
  {
    assert(m_initialized);
    for (unsigned int i = 0; i < size0 (); ++i)
    {
      m_data[i] *= a;
    }
    return *this;
  }

  template <typename Number>
  ForestVector<Number> &
  ForestVector<Number>::operator = (const Number a)
  {
    assert(m_initialized);
    for (unsigned int i = 0; i < size0 (); ++i)
    {
      m_data[i] = a;
    }
    return *this;
  }

  template <typename Number>
  void
  ForestVector<Number>::swap (ForestVector<Number> &v)
  {
    assert(m_initialized);
    Assert(m_n_multi_blocks == v.m_n_multi_blocks,
        ExcDimensionMismatch(m_n_multi_blocks, v.m_n_multi_blocks));
    for (unsigned int g = 0; g < m_n_multi_blocks; ++g)
    {
      m_data[g].swap (v.m_data[g]);
    }
  }

  template <typename Number>
  Number
  ForestVector<Number>::l2_norm () const
  {
    assert(m_initialized);
    Number val = 0;
    for (unsigned int i = 0; i < size0 (); ++i)
    {
      Number data_l2_norm = m_data[i].l2_norm ();
      val += data_l2_norm * data_l2_norm;
    }
    return std::sqrt (val);
  }

  template <typename Number>
  Number
  ForestVector<Number>::add_and_dot (
      const Number a,
      const ForestVector< Number > &V,
      const ForestVector< Number > &W)
  {
    assert(m_initialized);
    Number val = 0;
    for (unsigned int i = 0; i < size0 (); ++i)
    {
      val += m_data[i].add_and_dot (a, V.m_data[i], W.m_data[i]);
    }
    return val;
  }

  template <typename Number>
  std::size_t
  ForestVector<Number>::memory_consumption () const
  {
    std::size_t tmc = 0; /* total_memory_consumption */
    if (m_data.size () > 0)
    {
      tmc = m_data[0].memory_consumption () * m_data.size ();
    }
    return tmc;
  }

  template <typename Number>
  unsigned int
  ForestVector<Number>::size0 () const
  {
    return m_n_multi_blocks;
  }

  template <typename Number>
  const std::vector<unsigned int> &
  ForestVector<Number>::size1 () const
  {
    return m_block_sizes;
  }

  template class ForestVector<double> ;
  template class ForestVector<float> ;

//---------------------------------------------------------
/*

 template <class MATRIX>
 MyBlockMatrix<MATRIX>::MyBlockMatrix (unsigned int n_blocks)
 : block (n_blocks),
 n (n_blocks)
 {
 }

 template <class MATRIX>
 unsigned int
 MyBlockMatrix<MATRIX>::n_blocks ()
 {
 return n;
 }

 template class MyBlockMatrix<SparseMatrix<double> > ;

 template <class MATRIX, class VECTOR, class SVECTOR>
 OrigTransMatrix<MATRIX, VECTOR, SVECTOR>::OrigSource(unsigned int n)
 :
 nsigf(n),
 n(n)
 {}

 template <class MATRIX, class VECTOR, class SVECTOR>
 unsigned int OrigTransMatrix<MATRIX, VECTOR, SVECTOR>::get_n_blocks()
 {
 return n_blocks;
 }

 template <class MATRIX, class VECTOR, class SOURCE>
 void OrigTransMatrix<MATRIX, VECTOR, SVECTOR>::vmult (VECTOR &dst,
 const SOURCE &src) const
 {
 for (unsigned int i=0; i<n_blocks; ++i)
 {
 dst.block(i) = 0;
 for (unsigned int j=0; j<n_blocks; ++j)
 {
 scatt.block[j].vmult_add(dst.block(i), src.scalar);
 trans.block[i].vmult_add(dst.block(i), src.scalar);
 }
 }

 }

 template <class MATRIX, class VECTOR, class SVECTOR>
 void OrigTransMatrix<MATRIX, VECTOR, SVECTOR>::Tvmult (VECTOR &dst,
 const SOURCE &src) const
 {
 Assert (false, ExcNotImplemented());
 }

 template class OrigTransMatrix<SparseMatrix<double>,
 BlockVector<double>,
 OrigSource<SparseMatrix<double>,
 BlockVector<double>,
 Vector<double> > >;

 //---------------------------------------------------------

 template <class MATRIX, class VECTOR, class SOURCE>
 OrigTransMatrix<MATRIX, VECTOR, SOURCE>::OrigTransMatrix(unsigned int n)
 :
 trans(n),
 scatt(n),
 n(n)
 {}

 template <class MATRIX, class VECTOR, class SOURCE>
 unsigned int OrigTransMatrix<MATRIX, VECTOR, SOURCE>::get_n_blocks()
 {
 return n_blocks;
 }

 template <class MATRIX, class VECTOR, class SOURCE>
 void OrigTransMatrix<MATRIX, VECTOR, SOURCE>::vmult (VECTOR &dst,
 const SOURCE &src) const
 {
 for (unsigned int i=0; i<n_blocks; ++i)
 {
 dst.block(i) = 0;
 for (unsigned int j=0; j<n_blocks; ++j)
 {
 scatt.block[j].vmult_add(dst.block(i), src.scalar);
 trans.block[i].vmult_add(dst.block(i), src.scalar);
 } 
 }

 }

 template <class MATRIX, class VECTOR, class SOURCE>
 void OrigTransMatrix<MATRIX, VECTOR, SOURCE>::Tvmult (VECTOR &dst,
 const SOURCE &src) const
 {
 Assert (false, ExcNotImplemented());
 }

 template class OrigTransMatrix<SparseMatrix<double>,
 BlockVector<double>,
 OrigSource >;
 */

//---------------------------------------------------------
} // end of namespace Forest

