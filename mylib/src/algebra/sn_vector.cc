/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/sn_vector.cc
 * @brief  Implementation of classes SnVector and MyBlockMatrix
 */

#include "algebra/sn_vector.h"

#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/base/exceptions.h"    // for AssertDimension, Assert

namespace Forest
{
  using namespace dealii;

  SnVector::SnVector ()
      : ForestVector<double> ()
  {
  }

  SnVector::SnVector (const unsigned int n_groups,
                      const unsigned int n_ord,
                      const unsigned int n_elem)
      : ForestVector<double> (n_groups,
          std::vector<unsigned int> (n_ord, n_elem))
  {
  }

  SnVector::SnVector (const unsigned int n_groups,
                      const std::vector<unsigned int> & block_sizes)
      : ForestVector<double> (n_groups, block_sizes)
  {
  }

  void
  SnVector::reinit (const unsigned int n_groups,
                    const unsigned int n_ord,
                    const unsigned int n_elem,
                    bool leave_elements_uninitialized)
  {
    std::vector<unsigned int> block_sizes (n_ord, n_elem);
    BaseClass::reinit (n_groups, block_sizes, leave_elements_uninitialized);
  }


  /*void
  SnVector::reinit (const SnVector& v,
                    const bool leave_elements_uninitialized)
  {
    BaseClass::reinit (v, leave_elements_uninitialized);
  }*/

  /*double
  SnVector::operator * (const SnVector &v) const
  {
    assert(m_initialized);
    AssertDimension(size0(), v.size0());
    double val = BaseClass::operator* (v);
    return val;
  }*/

  /*
  void
  SnVector::add (const SnVector &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    BaseClass::add (x);
  }

  void
  SnVector::add (const double a,
                 const SnVector &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    BaseClass::add (a, x);
  }*/

  /*void
  SnVector::sadd (const double a,
                  const double b,
                  const SnVector &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    BaseClass::sadd (a, b, x);
  }*/

  /*
  void
  SnVector::equ (const double a,
                 const SnVector &x)
  {
    assert(m_initialized);
    AssertDimension(size0(), x.size0());
    BaseClass::equ (a, x);
  }
  */

} // end of namespace Forest

//---------------------------------------------------------

/* Create explicit instantiation for the vector class. If your project
 * consists of multiple files, including header files, this instantiation
 * must be put in a <code>.cc</code> file in order to instantiate only
 * once. */
#include "deal.II/lac/vector_memory.templates.h" // for GrowingVectorMemory, etc
template class dealii::VectorMemory<Forest::SnVector>;
template class dealii::GrowingVectorMemory<Forest::SnVector>;
