/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum.h
 * @brief  FissionSpectrum class template declarations
 */

#ifndef FOREST_FISSION_SPECTRUM_H
#define FOREST_FISSION_SPECTRUM_H

#include "algebra/sn_vector_fwd.h"

namespace Forest
{
  /**
   * @class FissionSpectrum
   * @ingroup ForestNeutronics
   * @brief Fission Spectrum operator \f$ \chi \f$
   *
   * @details
   * @remarks This class builds the matrix if requested, but also allows
   * matrix-free application of the fission source operator, both through
   * the same call to @ref vmult or @ref vmult_add
   * @todo update the documentation
   */
  template <int dim>
  class FissionSpectrum
  {

  public:

    /** @brief Virtual destructor */
    virtual
    ~FissionSpectrum ()
    {
    }
    ;

    /** @brief True when the operator is matrix free, false otherwise */
    virtual bool
    is_matrix_free () = 0;

    /**
     * @brief Apply the fission operator to @p phi and store
     * the result in @p source
     * @param source
     * @param phi
     */
    virtual void
    vmult (SnVector & source,
           const SnVector & phi) = 0;

    /**
     * @brief Apply the fission operator to @p phi and add
     * the result to @p source
     * @param source
     * @param phi
     */
    virtual void
    vmult_add (SnVector & source,
               const SnVector & phi) = 0;

    /**
     * @brief The memory required by the elements of the class
     * @return Memory consumption of object
     */
    virtual double
    memory_consumption () const = 0;

  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SPECTRUM_H */
