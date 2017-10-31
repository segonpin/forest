/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_spectrum_matrix_built.h
 * @brief  FissionSpectrumMatrixBuilt class template declarations
 */

#ifndef FOREST_FISSION_SPECTRUM_MATRIX_BUILT_H
#define FOREST_FISSION_SPECTRUM_MATRIX_BUILT_H

#include "neutronics/fission_spectrum.h"  // for FissionSpectrum

#include "neutronics/state_fwd.h"
#include "algebra/sn_vector_fwd.h"

#include "deal.II/lac/vector.h"         // for Vector

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim>  class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class FissionSpectrumMatrixBuilt
   * @brief Fission Spectrum operator \f$ \chi \f$, matrix-built version
   * @ingroup ForestNeutronics
   * @details @copydetails FissionSpectrum
   */
  template <int dim>
  class FissionSpectrumMatrixBuilt: public FissionSpectrum<dim>
  {

  public:

    /** @brief Constructor, receiving an associated state object */
    FissionSpectrumMatrixBuilt (State<dim> & state);

    /** @brief Destructor for this class. */
    virtual
    ~FissionSpectrumMatrixBuilt(){}

    /* Declaration of inherited pure virtual methods */

    virtual bool
    is_matrix_free ()
    {
      return false;
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

    void build();

    /**
     * @brief State of the problem containing the information needed.
     */
    State<dim> & mp_state;

    /**
     * @brief The fission spectrum matrix
     * @details This matrix is composed of diagonal matrices for each
     * group, when each matrix juts multiply by the value of \f$ \chi \f$
     * the dofs belonging to different cells for each energy.
     */
    std::vector<Vector<double> > m_chi_matrix;

    /**
     * @brief Reference to the DoFHandler object from @p state.
     */
    const DoFHandler<dim, dim> & mp_dof_handler;

    void set_non_zero_pattern();

    std::vector<bool> m_non_zero;

  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SPECTRUM_MATRIX_BUILT_H */
