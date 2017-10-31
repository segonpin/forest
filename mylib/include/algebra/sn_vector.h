/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/sn_vector.h
 * @brief   SnVector class template declarations
 */

#ifndef FOREST_SN_VECTOR_H
#define FOREST_SN_VECTOR_H

#include "algebra/forest_vector.h"

#include <vector>

/**
 * @brief  Classes related to algebra (matrix/vector)
 */

namespace Forest
{
  using namespace dealii;

  /**
   * @brief A vector class used in Forest
   * @ingroup ForestAlgebra
   *
   * @details This vector is used to hide some complexity.
   * This vector considers all the energy groups, and inside each energy
   * group it contains the scalar flux and the boundary conditions for each
   * direction. Then, we have to use it in different scenarios.
   *
   * First, we have to use it in a Krylov method as any other vector,
   * so the scaling functions, norms, and so on are defined for the whole
   * vector in a standard way way.
   *
   * Then we have to perform other operation only over the part corresponding
   * to the scalar flux, not to the boundary values. For this, we have added
   * the following constants:
   *
   * - n_lmom: number of Legendre moments to define the flux
   * - n_bv: number of angles for which there is an contribution from
   *   reflective boundary conditions.
   */
  class SnVector : public ForestVector<double>
  {
  public:
    typedef ForestVector<double> BaseClass;

    /**
     @brief Empty constructor.
     @details This is necessary when we want to create the object without
     specifying the dimensions. This should be initialize later.
     */
    SnVector ();

    /**
     @brief Standard constructor.
     @details Creates the object and allocate the necessary space.
     @param n_groups
     @param n_ord
     @param n_elem
     */
    SnVector (const unsigned int n_groups,
              const unsigned int n_ord,
              const unsigned int n_elem);

    /**
     @brief Standard constructor.
     @details Creates the object and allocate the necessary space.
     */
    SnVector (const unsigned int n_groups,
              const std::vector<unsigned int> & block_sizes);

  private:
    /** @brief False when called with empty constructor, true otherwise. */
    using BaseClass::m_initialized  ;

    /** @brief Number of groups in the current SnVector. */
    using BaseClass::m_n_multi_blocks;

    /** @brief Number of dofs each energy group. */
    using BaseClass::m_block_sizes;

  public:
    /** @brief data */
    using BaseClass::m_data;

    /**
     @brief Reinitialize the vector.
     @details Resize to have the structure provided by the arguments
     Note that the second argument must have
     a default value equal to false
     @param n_groups
     @param n_ord
     @param n_elem
     @param leave_elements_uninitialized = false
     */
    void
    reinit (const unsigned int n_groups,
            const unsigned int n_ord,
            const unsigned int n_elem,
            const bool leave_elements_uninitialized = false);

    /**
     @brief resize the vector.
     @details resize to have the same structure as the one provided and/or
     clear vector. note that the second argument must have
     a default value equal to false
     @param v
     @param leave_elements_uninitialized = false
    void
    reinit (const SnVector& v,
            const bool leave_elements_uninitialized = false);
     */

    using BaseClass::reinit;

    /**
     @brief Scalar product of this vector times vector @p v.
     @param v
     @retval res = <this,@p v>
    double
    operator * (const SnVector &v) const;
     */

    using BaseClass::operator*;

    /**
     @brief Addition of vectors.
     @details Add the vector @p x to the current vector.
     @code{.cpp} this += x @endcode
     @param x
    void
    add (const SnVector &x);
     */

    /**
     @brief Addition of scaled vector.
     @details Add the vector @p a @p x to the current vector.
     @code{.cpp} this += a*x @endcode
     @param a
     @param x
    void
    add (const double a,
         const SnVector &x);
     */

    using BaseClass::add;

    /**
     @brief Scaling and simple addition.
     @details Add the vector @p b @p x to the current vector scaled by @p a.
     @code{.cpp} this = a*this + b*x @endcode
     @param a
     @param b
     @param x
    void
    sadd (const double a,
          const double b,
          const SnVector &x);
     */

    using BaseClass::sadd;

    /**
     @brief Assignment.
     @details Assignment of a scaled vector.
     @code{.cpp} this = a*x @endcode
     @param a
     @param x
    void
    equ (const double a,
         const SnVector &x);
     */

    using BaseClass::equ;

    /**
     @brief Swap the content of this vector with the content of vector @p v.
     @details It only swaps the pointers to the data of the two vectors.
     @code{.cpp} tmp = this; this = v; v = tmp; @endcode
     @note They must have the same number of groups, directions and elements.
     @param v
    inline
    void
    swap (SnVector &v)
    {
      BaseClass::swap(v);
    }
     */

    using BaseClass::swap;

    using BaseClass::operator=;

    using BaseClass::memory_consumption;

    using BaseClass::size0;

    using BaseClass::size1;

    /** @brief Number of groups defining the first dimension. */
    inline
    unsigned int
    get_n_groups () const
    {
      return BaseClass::size0();
    }

    /** @brief Number of groups defining the first dimension. */
    inline
    const std::vector<unsigned int> &
    get_block_sizes () const
    {
      return BaseClass::size1();
    }

    // This is needed by gmres to use the right type for tolerances (float or double)
    typedef double value_type;

  };

} // end of namespace Forest

#endif /* FOREST_SN_VECTOR_H */
