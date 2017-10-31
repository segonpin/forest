/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/fission_source_matrix_built.h
 * @brief  FissionSourceMatrixBuilt class template declarations
 */

#ifndef FOREST_FISSION_SOURCE_MATRIX_BUILT_H
#define FOREST_FISSION_SOURCE_MATRIX_BUILT_H

#include "neutronics/fission_source.h"  // for FissionSource
#include "algebra/sn_vector_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#include "deal.II/lac/sparse_matrix.h"  // for SparseMatrix
#include "deal.II/lac/sparsity_pattern.h" // for SparsityPattern
#include "deal.II/meshworker/integration_info.h" // for IntegrationInfo
#include "deal.II/meshworker/dof_info.h"  // for DoFInfo

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim>  class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class FissionSourceMatrixBuilt
   * @ingroup ForestNeutronics
   * @brief Fission Source operator, matrix-built version
   *
   * @details @copydetails FissionSource
   *
   * @todo Add details about the matrix-built implementation
   */
  template <int dim>
  class FissionSourceMatrixBuilt: public FissionSource<dim>
  {

  public:

    /**
     * @brief Constructor, receiving an associated state object
     */
    FissionSourceMatrixBuilt (State<dim> & state);

    /** @brief Destructor for this class. */
    virtual
    ~FissionSourceMatrixBuilt()
    {
    }

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
     * @brief nsigf_integrate_cell_term
     * @param dinfo
     * @param info
     * @param g
     * @todo create a more general function depending on the parameter, and
     * inherit from this class for this matrix, the mass matrix,
     * the scattering matrix, and everything like this...
     */
    void
    nsigf_integrate_cell_term (MeshWorker::DoFInfo<dim> &dinfo,
                               MeshWorker::IntegrationInfo<dim> &info,
                               unsigned int g);

    /** @brief State of the problem containing the information needed. */
    State<dim> & mp_state;

    /** @brief Reference to the DoFHandler object from @p state. */
    const DoFHandler<dim, dim> & mp_dof_handler;

    SparsityPattern m_nusigf_pattern;

    /**
     * @brief The fission matrix.
     *
     * @details @p nsigf_matrix is composed of a product
     * of matrices, so in real we don't build the whole matrix
     * but the necessary for its application, as it is the different
     * fission matrices for the different energy groups and the
     * different fission spectra, to recover the whole matrix
     * in the following way
     *
     * \f{align*}{
     * \nu\Sigma_{f} = \left [
     * \begin{array}{ccc}
     * \ \chi_{0}I & \ldots & \chi_{G-1}I
     * \end{array}
     * \right ]
     * \cdot \left [
     * \begin{array}{c}
     *  \ \nu_{0}\Sigma_{f,0} \\
         *  \vdots \\
         *  \nu_{G-1}\Sigma_{f,G-1}
     * \end{array}
     * \right ]
     * \f}
     *
     */
    std::vector<SparseMatrix<double> > m_nusigf_matrix;

    void set_non_zero_pattern();

    std::vector<bool> m_non_zero;
  };

} // end of namespace Forest

#endif /* FOREST_FISSION_SOURCE_MATRIX_BUILT_H */
