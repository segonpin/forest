/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_bc.h
 * @brief  TransBc class template declarations
 */

#ifndef FOREST_TRANS_BC_H
#define FOREST_TRANS_BC_H

#include "algebra/sn_vector_fwd.h"
#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#ifndef DOXYGEN
namespace dealii { template <typename> class Vector; }
namespace dealii { template <typename> class BlockVector; }
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
#endif

#include <memory>

namespace Forest
{
  using namespace dealii;

  /**
   * @class TransBc
   * @ingroup ForestNeutronics
   * @brief Transport operator boundary contitions
   *
   * @details Defines the transport operator bountary conditions
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
   */

  template <int dim>
  class TransBc
  {
  public:

    /**
     * @brief Constructor for this class
     * @param state
     * @param ord
     */
    TransBc (State<dim> &state,
             std::shared_ptr<QuadratureBase<dim> > &ord);

    /**
     * @brief Add the reflective boundary conditions to the source term
     * @param source
     * @param phi_bc
     */
    void
    apply_LR (SnVector & source, const SnVector & phi_bc);

    /** @brief Add the reflective boundary conditions to the source term */
    void
    apply_LR (BlockVector<double> & source, const SnVector & phi_bc,
              unsigned int g);

    /** @brief Add the reflective boundary conditions to the source term */
    void
    apply_LR (Vector<double> & source, const SnVector & phi_bc, unsigned int g,
              unsigned int i_ord);

    /** @brief Add the reflective boundary conditions to the source term */
    void
    apply_LR (Vector<double> & source, const BlockVector<double> & phi_bc,
              unsigned int g, unsigned int i_ord);

    /**
     * @brief The same as @ref apply_LR but only for group g
     * @param source
     * @param sol_old
     * @param g
     */
    void
    apply_LR_g (BlockVector<double> & source,
                const BlockVector<double> & sol_old, unsigned int g);

    /**
     * @brief The same as @ref apply_LR_g but only for direction @p i_ord
     * @param source
     * @param sol_old
     * @param g
     * @param i_ord
     */
    void
    apply_LR_g_iord (Vector<double> & source,
                     const BlockVector<double> & sol_old, unsigned int g,
                     unsigned int i_ord);

    /** @todo Document me */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      return memory_consumption;
    }

  protected:

    /**
     * @brief The same as @ref apply_LR_g_iord for boundary @p boundary_face
     * @param dst
     * @param src
     * @param g
     * @param i_ord
     * @param boundary_face
     */
    void
    apply_LR_g_iord_f (Vector<double> & dst, const Vector<double> & src,
                       unsigned int g, unsigned int i_ord,
                       unsigned int boundary_face);

  private:

    const State<dim> & mp_state;
    const std::shared_ptr<QuadratureBase<dim> > & mp_ord;
    const DoFHandler<dim, dim> & mp_dof_handler;
    unsigned int m_n_angles;

  };

} // end of namespace Forest

#endif /* FOREST_TRANS_BC_H */
