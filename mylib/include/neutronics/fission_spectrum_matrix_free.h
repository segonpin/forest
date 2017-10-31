/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_matrix_free.h
 * @brief  FissionSpectrumMatrixFree class template declarations
 */

#ifndef FOREST_FISSION_SPECTRUM_MATRIX_FREE_H
#define FOREST_FISSION_SPECTRUM_MATRIX_FREE_H

#include "neutronics/fission_spectrum.h"  // for FissionSpectrum
#include "neutronics/state_fwd.h"
#include "algebra/sn_vector_fwd.h"

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim>  class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class FissionSpectrumMatrixFree
   * @ingroup ForestNeutronics
   * @brief Fission Spectrum operator \f$ \chi \f$, matrix-free version
   * @details @copydetails FissionSpectrum
   */
  template <int dim>
  class FissionSpectrumMatrixFree: public FissionSpectrum<dim>
  {

  public:

    /** @brief Constructor */
    FissionSpectrumMatrixFree (State<dim> & state);

    /** @brief Destructor for this class. */
    virtual
    ~FissionSpectrumMatrixFree(){}

    /* Declaration of inherited pure virtual methods */

    virtual bool
    is_matrix_free ()
    {
      return true;
    }

    virtual void
    vmult (SnVector & source,
           const SnVector & phi);

    virtual void
    vmult_add (SnVector & source,
               const SnVector & phi);

    virtual double
    memory_consumption () const;

  private:

    /**
     * @brief State of the problem containing the information needed.
     */
    State<dim> & mp_state;

    /**
     * @brief Reference to the DoFHandler object from @p state.
     */
    const DoFHandler<dim, dim> & mp_dof_handler;

    void set_non_zero_pattern();

    std::vector<bool> m_non_zero;

  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SPECTRUM_MATRIX_FREE_H */
