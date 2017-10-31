/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op.h
 * @brief  TransOp class template declarations
 */

#ifndef FOREST_TRANS_OP_H
#define FOREST_TRANS_OP_H

#include "algebra/sn_vector_fwd.h"

#ifndef DOXYGEN
namespace dealii { template <typename> class BlockVector; }
namespace dealii { template <typename> class Vector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class TransOp
   *
   * @ingroup ForestNeutronics
   *
   * @brief Transport operator and sweep, base class
   *
   * @details Defines the transport operator
   *
   * \f{align*}{
   * L(\mathbf{r},\Omega,E)\Psi(\mathbf{r},\Omega,E) =
   * \left( \Omega \cdot \nabla + \Sigma_{T}(\mathbf{r},E) \right)
   * \Psi(\mathbf{r},\Omega,E)
   * \f}
   *
   * with boundary conditions:
   *
   * - Incoming flux
   *
   *   \f{align*}{
   *   \Psi(\mathbf{r},\Omega,E) =
   *   \Psi_{in}(\mathbf{r},\Omega,E) \quad \forall \mathbf{r} \in \partial \Omega
   *   \f}
   *
   * - Albedo
   *
   *   \f{align*}{
   *   \Psi(\mathbf{r},\Omega,E) = \alpha
   *   \Psi_{in}(\mathbf{r},\Omega,E) \quad \forall \mathbf{r} \in \partial \Omega
   *   \f}
   *
   * and apply the transport sweep to the source vector \f$ Q(\mathbf{r},E) \f$
   *
   * \f{align*}{
   * D (L_0 + L_R) M \Phi = Q
   * \f}
   *
   * and this is apply performing several sweeps
   *
   * \f{align*}{
   * \Phi^{n+1} = D L_0^{-1} (M Q - L_R M \Phi^{n})
   * \f}
   *
   * \f{align*}{
   * \Phi^{n+1} = b - D L_0^{-1} L_R M \Phi^{n}
   * \f}
   *
   * where
   *
   * \f{align*}{
   * b = D L_0^{-1} M Q
   * \f}
   *
   * @note Different energy groups are completely uncoupled from the point
   * of view of the transport operator, so we can solve them separately.
   *
   * @todo update the documentation with weak formulation and matrix form
   */
  template <int dim>
  class TransOp
  {
  public:

    /** @brief Destructor */
    virtual
    ~TransOp ()
    {
    }

    /**
     * @brief Apply the inverse of the transport operator (the transport sweep)
     *
     * @param phi
     * @param source
     */
    void
    apply_L0_inv (SnVector & phi,
                  const SnVector & source);

    /** @brief Apply the inverse of the transport operator (the transport sweep) */
    void
    apply_L0_inv (SnVector & phi,
                  const BlockVector<double> & source,
                  const unsigned int g);

    /** @brief Apply the inverse of the transport operator (the transport sweep) */
    void
    apply_L0_inv (SnVector & phi,
                  const Vector<double> & source,
                  const unsigned int g,
                  const unsigned int i_ord);

    /** @brief Apply the inverse of the transport operator (the transport sweep) */
    void
    apply_L0_inv (Vector<double> & phi,
                  const Vector<double> & source,
                  const unsigned int g,
                  const unsigned int i_ord)
    {
      apply_L0_inv_g_iord (phi, source, g, i_ord);
    }

    /** @brief Apply the inverse of the transport operator (the transport sweep) */
    void
    apply_L0_inv_g (BlockVector<double> & phi,
                    const BlockVector<double> & source,
                    const unsigned int g);

    // pure virtual functions

    /** @brief Apply the inverse of the transport operator (the transport sweep) */
    virtual
    void
    apply_L0_inv_g_iord (Vector<double> & phi,
                         const Vector<double> & source,
                         const unsigned int g,
                         const unsigned int i_ord) = 0;

    /** @brief Access function */
    virtual
    unsigned int
    get_n_groups () const = 0;

    /** @brief Access function */
    virtual
    unsigned int
    get_n_angles () const = 0;

    /** @brief Returns object memory consumption */
    virtual
    double
    memory_consumption () const = 0;

    /** @brief access function */
    virtual
    double
    get_sweep_cputime() const = 0;

    virtual
    double
    get_sweep_walltime() const = 0;
  };


} // end of namespace Forest

#endif /* FOREST_TRANS_OP_H */
