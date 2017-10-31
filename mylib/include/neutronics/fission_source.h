/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source.h
 * @brief  FissionSource class template declarations
 */

#ifndef FOREST_FISSION_SOURCE_H
#define FOREST_FISSION_SOURCE_H

#include "algebra/sn_vector_fwd.h"

namespace Forest
{
  /**
   * @class FissionSource
   * @ingroup ForestNeutronics
   * @brief Fission Source operator, base class
   *
   * @details Defines the isotropic fission operator to be applied to a flux vector
   *
   * \f{align*}{
   * Q(\mathbf{r},E,t) = \frac{\chi(E)}{4\pi}\int_{0}^{\infty}\nu(E^\prime)\Sigma_{f}(\mathbf{r},E^{\prime},t)
   * \Phi(\mathbf{r},E^{\prime},t) dE^{\prime}
   * \f}
   *
   * We prepare different functions to use the fission source operator in a matrix-free fashion,
   * allowing to work with a single energy group
   *
   * \f{align*}{
   * Q_{g}(\mathbf{r},t) = \frac{\chi_{g}}{4\pi}\sum_{h \in G}\nu_{h}\Sigma_{f,h}(\mathbf{r},t)
   * \Phi_{h}(\mathbf{r},t)
   * \f}
   *
   * @remarks This class builds the matrix if requested, but also allows
   * matrix-free application of the fission source operator, both through
   * the same call to @ref vmult or @ref vmult_add
   *
   * @todo update the documentation
   *
   */
  template <int dim>
  class FissionSource
  {

  public:

    /** @brief Virtual destructor for this class. */
    virtual
    ~FissionSource()
    {
    };

    /** @brief True when the operator is matrix free, false otherwise. */
    virtual bool is_matrix_free() = 0;

    /**
     * @brief Apply the fission operator to @p phi and store
     * the result in @p source
     *
     * @param source
     * @param phi
     */
    virtual void
    vmult (SnVector & source,
           const SnVector & phi) = 0;

    /**
     * @brief Apply the fission operator to @p phi and add
     * the result to @p source
     *
     * @param source
     * @param phi
     */
    virtual void
    vmult_add (SnVector & source,
               const SnVector & phi) = 0;

    /**
     * @brief The memory required by the elements of the class.
     *
     * @return Memory consumption of object
     */
    virtual double
    memory_consumption () const = 0;

  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SOURCE_H */
