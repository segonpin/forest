/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op.h
 * @brief  DiffOp class template declarations
 */

#ifndef FOREST_DIFF_OP_H
#define FOREST_DIFF_OP_H

#include "algebra/sn_vector_fwd.h"

#ifndef DOXYGEN
namespace dealii { template <typename > class BlockVector; }
namespace dealii { template <typename > class Vector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class DiffOp
   * @ingroup ForestNeutronics
   * @brief Diffusion operator and sweep, base class
   *
   * @details Defines the diffusion operator
   *
   *   \f{align*}{
   *   & - \nabla \frac{1}{3\Sigma_{T,g}(\mathbf{r})} \nabla \Phi_g(\mathbf{r}) 
   *    + \Sigma_{T,g}(\mathbf{r}) \Phi_g(\mathbf{r})  \nonumber \   \
   *    & = \frac{1}{4\pi} \sum_{h=0}^{G-1}
   *    \Sigma_{s0,h\to g}(\mathbf{r}) \Phi_{h}(\mathbf{r})
   *    + \frac{1}{\lambda}\chi_g(\mathbf{r})\frac{1}{4\pi}\sum_{h=0}^{G-1}
   *    \nu\Sigma_{f,h}(\mathbf{r})
   *    \Phi_h(\mathbf{r}) 
   *    \f}
   *
   * or using the diffusion coefficient and the absorption cross section
   *
   *    \f{align*}{
   *    & - \nabla D_{g}(\mathbf{r}) \nabla \Phi_g(\mathbf{r}) 
   *    + \Sigma_{a,g}(\mathbf{r}) \Phi_g(\mathbf{r})  \nonumber \   \
   *    & = \frac{1}{4\pi} \sum_{h=0,h\neq g}^{G-1}
   *    \Sigma_{s0,h\to g}(\mathbf{r}) \Phi_{h}(\mathbf{r})
   *    + \frac{1}{\lambda}\chi_g(\mathbf{r})\frac{1}{4\pi}\sum_{h=0}^{G-1}
   *    \nu\Sigma_{f,h}(\mathbf{r})
   *    \Phi_h(\mathbf{r}) 
   *    \f}
   *
   * @todo Correct the documentation
   */
  template <int dim>
  class DiffOp
  {
  public:

    /** @brief Destructor */
    virtual
    ~DiffOp ()
    {
    }

    /**
     * @brief Apply the diffusion operator - SnVector version
     * @param phi
     * @param source
     */
    void
    vmult (SnVector & phi,
           const SnVector & source);

    /**
     * @brief Apply the diffusion operator for one group - SnVector version
     * @param phi
     * @param source
     * @param g
     */
    void
    vmult (SnVector & phi,
           const BlockVector<double> & source,
           const unsigned int g);

    // pure virtual functions

    /**
     * @brief Apply the diffusion operator for one group - BlockVector version
     * @param phi
     * @param source
     * @param g
     */
    virtual
    void
    vmult (BlockVector<double> & phi,
           const BlockVector<double> & source,
           const unsigned int g) = 0;

    /** @brief Preconditioner ? */
    virtual
    void
    prec (BlockVector<double> & dst,
          const BlockVector<double> & src,
          const unsigned int g) const = 0;

    /** @brief Access function */
    virtual
    unsigned int
    get_n_groups () const = 0;

    /** @brief Access function */
    virtual
    double
    memory_consumption () const = 0;
  };

} // end of namespace Forest

#endif /* FOREST_DIFF_OP_H */
