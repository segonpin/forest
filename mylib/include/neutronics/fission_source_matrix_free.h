/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_matrix_free.h
 * @brief  FissionSourceMatrixFree class template declarations
 */

#ifndef FOREST_FISSION_SOURCE_MATRIX_FREE_H
#define FOREST_FISSION_SOURCE_MATRIX_FREE_H

#include "neutronics/fission_source.h"  // for FissionSource
#include "algebra/sn_vector_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim>  class DoFHandler; }
namespace dealii { template <typename Number> class Vector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class FissionSourceMatrixFree
   * @ingroup ForestNeutronics
   * @brief Fission Source operator, matrix-free version
   *
   * @details @copydetails FissionSource
   *
   * @todo Add details about the matrix-free implementation
   */
  template <int dim>
  class FissionSourceMatrixFree: public FissionSource<dim>
  {

  public:

    /**
     * @brief Constructor, receiving an associated state object
     */
    FissionSourceMatrixFree (State<dim> & state);


    /** @brief Destructor */
    virtual
    ~FissionSourceMatrixFree()
    {
    }

    /* Declaration of inherited pure virtual methods */

    virtual bool
    is_matrix_free ()
    {
      return true;
    }

    void
    vmult (SnVector & source,
           const SnVector & phi);

    void
    vmult_add (SnVector & source,
               const SnVector & phi);

    double
    memory_consumption () const;

  private:

    /**
     * @brief calculate the fission source for group @p g
     * and add it to @p source
     * @details The local application of the finite  element operator
     *
     * \f{align*}{
     * \mathcal{M} (u,v) = \int_{\Omega_e} c_{e}
     * u(\mathbf{r})v(\mathbf{r})
     * \f}
     *
     * where
     *
     * \f{align*}{
     * c_{e} \nu_{g}\Sigma_{f,g}(\mathbf{r})
     * \f}
     *
     * @param source
     * @param phi
     * @param g
     *
     * @todo Prepare sum-factorization technique
     */
    void
    ingroup_fission_add (Vector<double> & source,
                         const Vector<double> & phi,
                         unsigned int g);

    /** @brief State of the problem containing the information needed. */
    State<dim> & mp_state;

    /** @brief Reference to the DoFHandler object from @p state. */
    const DoFHandler<dim, dim> & mp_dof_handler;

    void set_non_zero_pattern();

    std::vector<bool> m_non_zero;

  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SOURCE_MATRIX_FREE_H */
