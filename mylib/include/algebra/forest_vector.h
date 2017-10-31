/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/forest_vector.h
 * @brief   ForestVector (and MyBlockMatrix) class template declaration(s)
 */

#ifndef FOREST_VECTOR_H
#define FOREST_VECTOR_H

#include "utils/forest_vector.templates.h"
#include "deal.II/lac/block_vector.h"

#include <vector>
#include <cstddef>

/**
 * @brief  Classes related to algebra (matrix/vector)
 */
namespace Forest
{
  using namespace dealii;

  /**
   * @class ForestVector
   * @ingroup ForestAlgebra
   * @brief A vector class used in Forest
   */
  template <typename Number>
  class ForestVector
  {
  public:

    /**
     @brief Empty constructor.
     @details This is necessary when we want to create the object without
     Specifying the dimensions. This should be initialize later.
     */
    ForestVector ();

    /**
     @brief Standard constructor.
     @details Creates the object and allocate the necessary space.
     */
    ForestVector (const unsigned int n_multi_blocks,
                  const std::vector<unsigned int> & block_sizes);

  protected:

    /** @brief False when called with empty constructor, true otherwise. */
    bool m_initialized;

    /** @brief Number of multi-blocks in the current Vector. */
    unsigned int m_n_multi_blocks;

    /** @brief Number of dofs each energy group. */
    std::vector<unsigned int> m_block_sizes;

    /** @brief data */
    std::vector<BlockVector<Number> > m_data;

  public:

    /**
     @brief Reinitialize the vector.
     @details Resize to have the structure provided by the arguments
     Note that the second argument must have
     a default value equal to false
     @param n_multi_blocks
     @param block_sizes
     @param leave_elements_uninitialized = false
     */
    void
    reinit (const unsigned int n_multi_blocks,
            const std::vector<unsigned int> & block_sizes,
            const bool leave_elements_uninitialized = false);

    /**
     @brief resize the vector.
     @details resize to have the same structure as the one provided and/or
     clear vector. note that the second argument must have
     a default value equal to false
     @param v
     @param leave_elements_uninitialized = false
     */
    void
    reinit (const ForestVector<Number>& v,
            bool leave_elements_uninitialized = false);

    /**
     @brief Scalar product of this vector times vector @p v.
     @param v
     @retval res = <this,@p v>
     */
    Number
    operator * (const ForestVector<Number> &v) const;

    /**
     @brief Addition of vectors.
     @details Add the vector @p x to the current vector.
     @code{.cpp} this += x @endcode
     @param x
     */
    void
    add (const ForestVector<Number> &x);

    /**
     @brief Addition of scaled vector.
     @details Add the vector @p a @p x to the current vector.
     @code{.cpp} this += a*x @endcode
     @param a
     @param x
     */
    void
    add (const Number a,
         const ForestVector<Number> &x);

    /**
     @brief Scaling and simple addition.
     @details Add the vector @p b @p x to the current vector scaled by @p a.
     @code{.cpp} this = a*this + b*x @endcode
     @param a
     @param b
     @param x
     */
    void
    sadd (const Number a,
          const Number b,
          const ForestVector<Number> &x);

    /**
     @brief Assignment.
     @details Assignment of a scaled vector.
     @code{.cpp} this = a*x @endcode
     @param a
     @param x
     */
    void
    equ (const Number a,
         const ForestVector<Number> &x);

    /**
     @brief Scale the elements of the vector by a fixed value.
     @details Scaling this vector with the scalar value @p a.
     @code{.cpp} this *= a @endcode
     @param a
     */
    ForestVector<Number> &
    operator *= (const Number a);

    /**
     @brief Assignment.
     @details Assign the vector the constant value @p a.
     @code{.cpp} this = a @endcode
     @param a
     */
    ForestVector<Number> &
    operator = (const Number a);

    /**
     @brief Swap the content of this vector with the content of vector @p v.
     @details It only swaps the pointers to the data of the two vectors.
     @code{.cpp} tmp = this; this = v; v = tmp; @endcode
     @note They must have the same number of groups, directions and elements.
     @param v
     */
    void
    swap (ForestVector<Number> &v);

    /**
     @brief Return the l2 norm of the vector.
     @retval res = <this,this>^{1/2}
     */
    Number
    l2_norm () const;

    /**
     Performs a combined operation of a vector addition and a
     subsequent inner product, returning the value of the inner product.
     */
    Number
    add_and_dot (const Number a, const ForestVector< Number > &V, const ForestVector< Number > &W);

    /** @brief memory_consumption. */
    std::size_t
    memory_consumption () const;

    /** @brief Number of groups defining the first dimension. */
    unsigned int
    size0 () const;

    /** @brief Number of directions defining the second dimension. */
    const std::vector<unsigned int> &
    size1 () const;
  };

/*
 template <class MATRIX>
 class MyBlockMatrix
 {
 public:
 MyBlockMatrix (unsigned int n);

 unsigned int
 n_blocks ();

 std::vector<MATRIX> block;

 private:
 unsigned int n;
 };

 template <class MATRIX, class VECTOR, class SVECTOR>
 class OrigSource
 {
 public:
 /// Constructor
 OrigSource(unsigned int n);
 /// Application of matrix to vector src.
 /// Write result into dst
 void init (const VECTOR &dst) const;
 /// Application of transpose to a Vector.
 /// Only used by certain iterative methods.
 void update(unsigned int i_block, const VECTOR &dst) const;
 unsigned int n_blocks();
 MyBlockMatrix<MATRIX> &nsigf;
 VECTOR src_i;
 SVECTOR src;
 private:
 unsigned int n;
 };
 template <class MATRIX, class VECTOR, class SOURCE>
 class OrigTransMatrix
 {
 public:
 /// Constructor
 OrigTransMatrix(unsigned int n);
 /// Application of matrix to vector src.
 /// Write result into dst
 void vmult (VECTOR &dst, const SOURCE &src) const;
 /// Application of transpose to a Vector.
 /// Only used by certain iterative methods.
 void Tvmult (VECTOR &dst, const SOURCE &src) const;
 unsigned int n_blocks();
 MyBlockMatrix<MATRIX> &trans;
 MyBlockMatrix<MATRIX> &scatt;
 private:
 unsigned int n;
 };
 */

} // end of namespace Forest

#endif /* FOREST_VECTOR_H */
