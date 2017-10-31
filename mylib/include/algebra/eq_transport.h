/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eq_transport.h
 * @brief   EqTransport class template declarations
 */

#ifndef FOREST_EQ_TRANSPORT_H
#define FOREST_EQ_TRANSPORT_H

#include "algebra/equation.h"             // for SnVector, Equation
#include "algebra/solver_eq.h"            // for EqSolverWG and EqPrec
#include "angle/directions_to_moments.h"  // for DirectionsToMoments
#include "angle/moments_to_directions.h"  // for MomentsToDirections
#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"         // For class State&

#include "deal.II/base/exceptions.h"      // for Assert
#include "deal.II/lac/block_vector.h"     // for BlockVector
#include "deal.II/lac/vector.h"           // for Vector

#include <memory>                         // for shared_ptr
#include <vector>                         // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class ScatteringSource; }
namespace Forest { template <int dim> class FissionSource; }
namespace Forest { template <int dim> class FissionSpectrum; }
namespace Forest { template <int dim> class BoundaryValues; }
namespace Forest { template <int dim> class TransOp; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class EqTransport
   * @ingroup ForestAlgebra
   * @brief Neutron transport equation
   *
   * @details This class implements the transport equation at two
   * different levels. We have to solve a source problem in a multi-group
   * level, and this is approximated by working with single-group solvers.
   *
   * # Multi-Group (Gauss-Seidel)
   *
   * We have to solve the following source problem
   *
   *  \f{align*}{
   *  L \Psi = MSD \Psi + Q
   *  \f}
   *
   * Where the source \f$ Q \f$ may contain the fission source scaled with the
   * multiplication factor and other external sources. This system has a block
   * structure defined by the energy discretization, where the streaming
   * operator \f$ L \f$ is block diagonal (no coupling between different
   * groups), and operator \f$ S \f$ is "almost" lower triangular
   * (up-scattering is usually smaller than down-scattering). A reasonable
   * approach would be then to solve the systems as if it would be lower
   * triangular, and "lagger" the up-scattering to be multiplied to a guess.
   * To do so, we start solving the system for the first energy group, and then
   * applying the Gauss-Seidel method as follows
   *
   * \f{align*}{
   * (L_g+MS_{g \to g}D)\Psi_{g}^{n+1} =
   * \sum_{g'=1}^{g-1}M S_{g' \to g} D \Psi_{g'}^{n+1} +
   * \sum_{g'=g+1}^{G}M S_{g' \to g} D \Psi_{g'}^{n} + Q_{g}
   * \quad g = 1, \ldots, G .
   * \f}
   *
   * In order to simplify the notation, grouping right hand side in a unique
   * source per group
   *
   * \f{align*}{
   * & \hat{Q}_{g} = \sum_{g'=1}^{g-1}M S_{g' \to g} D \Psi_{g'}^{n+1} +
   * \sum_{g'=g+1}^{G}M S_{g' \to g} D \Psi_{g'}^{n} + Q_{g}
   * \quad g = 1, \ldots, G .
   * \f}
   *
   * we rewrite the system as
   *
   * \f{align*}{
   * & (L_g+MS_{g \to g}D)\Psi_{g}^{n+1} = \hat{Q}_{g}
   * \quad g = 1, \ldots, G
   * \f}
   *
   * @note It is important to notice here that using this strategy we are
   * not really solving the system. In order to do so, this iteration can be
   * used as a preconditioner for a full energy groups GMRES. The
   * update obtained here is a good approximation, and it is enough when using
   * the power iteration to solve the eigenvalue problem, while it would not
   * be enough if the eigenproblem is solved with, for instance, an Arnoldi
   * method.
   *
   * # Single Group Source Iteration.
   *
   * Here we perform, for a given energy group @p m_g, one
   * source iteration. Lets consider the following problem (we skip the group
   * index \f$ g \f$ because it happens inside each group)
   *
   * \f{align*}{
   * L \Psi = MSD \Psi + Q
   * \f}
   *
   * Of course we can directly try to invert the matrix \f$ (L-MSD)\f$ and get
   * the solution \f$ \Psi = (L-MSD)^{-1} Q\f$, but this system is very large,
   * and its structure suggests us a more efficient procedure, called source
   * iteration.
   *
   * Now, the streaming operator can be decomposed into the part without
   * reflecting boundary conditions, \f$ L_{0} \f$ , and the coupling between
   * the different directions due to the reflecting boundary conditions,
   * \f$ L_{R} \f$ , as follows
   *
   * \f{align*}{
   * L = L_{0} + L_{R}
   * \f}
   *
   * Then, two scenarios are considered. The first one neglect reflective
   * boundary conditions, in order to simplify the explanation, and the second
   * one includes reflecting boundary conditions.
   *
   * ## Neglecting Reflecting Boundary Conditions
   *
   * Here we assume \f$ L_{R} \equiv 0 \f$, so \f$ L = L_{0} \f$ is a block
   * diagonal matrix that can be easily inverted by performing transport
   * sweeps. We "lagged" the scattering term to a previous iteration step,
   * we use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality, and then update the solution inverting the
   * streaming term (what should be easier)
   *
   * \f{align*}{
   * & \Psi^{l+1} = L^{-1}(MS \Phi^{l} + Q) \\
   * & \Phi^{l} = D \Psi^{l}
   * \f}
   *
   * A useful trick here, is to use, instead of the source \f$ Q  \f$, what is
   * called the uncollided source \f$ b = L^{-1} Q \f$, and then perform the
   * following iteration procedure
   *
   * \f{align*}{
   * & \Psi^{l+1} = L^{-1} MS \Phi^{l} + b \\
   * & \Phi^{l} = D \Psi^{l}
   * \f}
   *
   * This expression, ignoring boundary conditions, can be simplified further
   * by combining the two equations as follows
   *
   * \f{align*}{
   * & \Phi^{l+1} = D L^{-1} MS \Phi^{l} + D b
   * \f}
   *
   * what is equivalent to solve the system
   *
   * \f{align*}{
   * & (I - D L^{-1} MS) \Phi = D b
   * \f}
   *
   *
   * ## Including Reflecting Boundary Conditions
   *
   *
   * Here we assume "lagged" the scattering term to a previous iteration step,
   * as before, but also "lagged" the effect of the boundary conditions, and
   * we use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality. Thus we update the solution inverting the
   * streaming term (block diagonal)
   *
   * \f{align*}{
   * & \Psi^{l+1} = L_{0}^{-1}(MS \Phi^{l} - L_{R}\Psi^{l}) + b \\
   * & \Phi^{l} = D \Psi^{l}
   * \f}
   *
   * where here \f$ b \f$ is defined as the uncollided flux
   * \f$ b := L_{0}^{-1} Q \f$, without considering the boundary conditions.
   * Now we can not remove the angular flux combining the two equations
   * because of the boundary conditions. Instead, we are going to consider the
   * boundary conditions applying to the vector of significant dofs, that is
   * projected to the full vector through the operator \f$ P \f$ , and we
   * consider our state vector the vector of the scalar flux and the vector
   * of significant angular dofs, as follows
   *
   * \f{align*}{
   * & \Psi^{l+1} = L_{0}^{-1}(MS \Phi^{l} - L_{R}P\Psi_{S}^{l}) + b \\
   * & \Phi^{l} =  D \Psi^{l} \\
   * & \Psi_{S}^{l+1} = P^{T}\Psi^{l+1}
   * \f}
   *
   * that can be reduced to
   *
   * \f{align*}{
   * & \Phi^{l+1} = D L_{0}^{-1}(MS \Phi^{l} - L_{R}P\Psi_{S}^{l}) + D b \\
   * & \Psi_{S}^{l+1} = P^{T}L_{0}^{-1}(MS \Phi^{l} - L_{R}P\Psi_{S}^{l}) + P^{T}b
   * \f}
   *
   * Or in a matrix form, as follows
   *
   * \f{align*}{
   * \left[ \begin{array}{c}
   * \Phi^{l+1} \\
   * \Psi^{l+1}_{S}
   * \end{array} \right]
   * =
   * \left[ \begin{array}{cc}
   * D & 0 \\
   * 0 & I
   * \end{array} \right]
   * \left(
   * \left[ \begin{array}{c}
   * I \\
   * P^{T}
   * \end{array} \right]
   * L_{0}^{-1}
   * \left[ \begin{array}{cc}
   * I & -L_{R} P
   * \end{array} \right]
   * \right)
   * \left[ \begin{array}{cc}
   * M & 0 \\
   * 0 & I
   * \end{array} \right]
   * \left[ \begin{array}{cc}
   * S & 0 \\
   * 0 & I
   * \end{array} \right]
   * \left[ \begin{array}{c}
   * \Phi^{l} \\
   * \Psi^{l}_{S}
   * \end{array} \right]
   * +
   * \left[ \begin{array}{c}
   * D b \\
   * P^{T}b
   * \end{array} \right]
   * \f}
   *
   * This can be rewritten as
   *
   * \f{align*}{
   * & \hat{\Phi}^{l+1} = \hat{D} \hat{L}^{-1} \hat{M}\hat{S} \hat{\Phi}^{l}
   * + \hat{b}
   * \f}
   *
   * where
   *
   * \f{align*}{
   * \hat{\Phi} :=
   * \left[ \begin{array}{c}
   * \Phi\\
   * \Psi_{S}
   * \end{array} \right],
   * \quad \hat{L}^{-1} :=
   * \left(
   * \left[ \begin{array}{c}
   * I \\
   * P^{T}
   * \end{array} \right]
   * L_{0}^{-1}
   * \left[ \begin{array}{cc}
   * I & -L_{R} P
   * \end{array} \right]
   * \right)
   * \equiv
   * \left[ \begin{array}{cc}
   * L_{0}^{-1}      & -L_{0}^{-1} L_{R} P \\
   * P^{T}L_{0}^{-1} & -P^{T}L_{0}^{-1} L_{R} P
   * \end{array} \right]
   * \f}
   *
   * and the rest of the operators follow the notation in the obvious way.
   *
   * @attention In the paper @cite warsa2004krylov the vector \f$ \hat{b} \f$ has
   * the second component equal to zero, while here I obtain that the second
   * component is equal to \f$ P^{T} b \f$. Check if it should be the case.
   *
   * We can also rearrange this system to use a Krylov solver as follows
   *
   * \f{align*}{
   * & (\hat{I} - \hat{D} \hat{L}^{-1} \hat{M}\hat{S}) \hat{\Phi} = \hat{b}
   * \f}
   *
   */
  template <int dim>
  class EqTransport : public Equation<dim>
  {
  public:

    /**
     * @brief Constructor of the Neutron Transport Equation.
     * @details It initializes all the operators (L0, LR, S, F, M, D).
     * @param state
     */
    EqTransport (State<dim> & state);

    /** @brief Destructor. */
    virtual
    ~EqTransport ()
    {
    }

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_angular_sizes () const;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_flux_bv_sizes () const;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_fission_density_sizes () const;

    /** @brief Is the system symmetric?. */
    virtual bool
    is_symmetric ()
    {
      return false;
    }


    /** @todo document.  */
    virtual void
    vmult (SnVector & dst,
           const SnVector & src) const;

    /** @brief Preconditioner ? */
    virtual void
    prec (BlockVector<double> & dst,
          const BlockVector<double> & src) const;

    /** @brief set the Preconditioner */
    virtual void
    set_prec (std::shared_ptr<Equation<dim> > & prec) const;


    virtual void
    solve (BlockVector<double> &dst,
           const BlockVector<double> &src,
           const double tol_default = -1.,
           const unsigned int max_it_default = 0) const;

    /**
     * @brief Application of matrix to vector src, and write result into dst.
     * @details Here we perform, for a given energy group @p m_g, one
     * source iteration. In a compact form it can be expressed as follows
     * \f{align*}{
     * & \hat{x} = (\hat{I} - \hat{D} \hat{L}^{-1} \hat{M}\hat{S}) \hat{b}
     * \f}
     * @param dst
     * @param src
     */
    virtual void
    vmult (BlockVector<double> &dst,
           const BlockVector<double> &src) const;

    /** @brief access to the scattering method. */
    /*virtual void
    scatt_add (BlockVector<double> &dst,
               const BlockVector<double> &src,
               const unsigned int g,
               const unsigned int h) const;*/

    /** @brief access to the up-scattering method. */
    virtual void
    up_add (BlockVector<double> &dst,
            const SnVector &src,
            const unsigned int g) const;

    /** @brief access to the down-scattering method. */
    virtual void
    down_add (BlockVector<double> &dst,
              const SnVector &src,
              const unsigned int g) const;

    /** @brief access to the fission method. */
    /* virtual void
    update_fission_source (SnVector &dst,
                           const SnVector &src); */

    /** @brief access to the fission method. */
    virtual void
    calculate_fission_density (SnVector &dst,
                               const SnVector &src);

    /** @brief access to the fission method. */
    virtual void
    distribute_fission_rate (SnVector &dst,
                             const SnVector &src);

    /**
     * @brief Application of the transpose of a matrix to a vector.
     * Used by some iterative methods.
     */
    virtual void
    Tvmult (SnVector & /* dst */,
            const SnVector & /* src */) const
    {
      Assert(false, ExcNotImplemented());
    }

    /**
     * @brief Application of the transpose of a matrix to a vector.
     * Used by some iterative methods.
     */
    virtual void
    Tvmult (BlockVector<double> & /* dst */,
            const BlockVector<double> & /* src */) const
    {
      Assert(false, ExcNotImplemented());
    }

    /**
     * @brief Preparing the right hand side for the source iteration.
     * @details For a given energy group @p m_g, we generate the right
     * hand side as follows
     *
     * \f{align*}{
     * \left[ \begin{array}{c}
     * x_{0} \\
     *  x_{S}
     * \end{array} \right]
     * =
     * \left[ \begin{array}{cc}
     * D & 0 \\
     *  0 & I
     * \end{array} \right]
     * \left[ \begin{array}{c}
     * I \\
     *  P^{T}
     * \end{array} \right]
     * L_{0}^{-1} M Q
     * \f}
     *
     * @param dst
     * @param src
     */
    virtual void
    set_rhs (BlockVector<double> & dst,
             const BlockVector<double> &src) const;

    virtual void
    get_angular (SnVector & dst,
                 const SnVector & src,
                 const double keff) const;

    /** @brief Set the group we are going to perform the group iteration to. */
    virtual void
    set_group (const unsigned int g)
    {
      m_g = g;
    }

    virtual double
    memory_consumption () const;

    /** @brief Set the group we are going to perform the group iteration to. */
    virtual unsigned int
    get_n_angles () const;

    virtual unsigned int
    get_n_groups() const;

    /** @brief Return the total number of sweeps within energy groups. */
    unsigned int
    get_tot_wg_it () const
    {
      return m_solver.get_tot_wg_it ();
    }

  private:

    /**
     * @brief Application of matrix to vector src, and write result into dst.
     * @details Here we perform, for a given energy group @p m_g, one
     * source iteration. In a matrix form it can be expressed as follows
     *
     * \f{align*}{
     * & \hat{x} = \hat{D} \hat{L}^{-1} \hat{M}\hat{S} \hat{b}
     * \f}
     *
     * We use the following notation for the arguments:
     *
     * \f{align*}
     * \text{src} := b \equiv
     * \left[ \begin{array}{c} b_{0} \\ b_{S} \end{array} \right]
     * \in \mathbb{R}^{I \times M} \times \mathbb{R}^{S} ,
     * \quad \text{and} \quad
     * \text{dst} := x \equiv
     * \left[ \begin{array}{c} x_{0} \\ x_{S} \end{array} \right]
     * \in \mathbb{R}^{I \times M} \times \mathbb{R}^{S}
     * \f}
     *
     * Where \f$ I \f$ is the number of degrees of freedom for the spatial
     * discretization over the mesh, \f$ M \f$ is the number of moments in
     * the scattering operator, and \f$ S \f$ is the number of significant
     * degrees of freedom for the boundary conditions.
     *
     * We perform the following steps:
     *
     * - The first thing is to initialize the vector \f$ x \f$ to zero
     * \f{align*}
     * x :=
     * \left[ \begin{array}{c} x_{0} \\ x_{S} \end{array} \right] =
     * \left[ \begin{array}{c} 0     \\ 0 \end{array} \right]
     * \f}
     *
     * - Then we put the scattering generated into m_rhs
     * \f{align*}
     * v_{1} = S b_{0}
     * \f}
     *
     * - Now, for every direction @p n we do the following operations
     * \f{align*}
     * & \text{for n = 1:N} &  \\
     * &&& v_{2} = M_{n} v_{1} \\
     * &&& v_{2} = v_{2} - P L_{R,n} b_{S,n} \\
     * &&& v_{3} = L_0^{-1} v_{2} \\
     * &&& x_{S} = x_{S} + P^{T} v_{3} \\
     * &&& x_{0} = x_{0} + D_{n} v_{3} \\
     * & \text{endfor} &
     * \f}

     * @param dst
     * @param src
     */
    void
    DL_invMS (BlockVector<double> &dst,
              const BlockVector<double> &src) const;

    void
    DL_invMS (SnVector &dst,
              const SnVector &src) const;

    void
    DL_invM (BlockVector<double> &dst,
             const BlockVector<double> &src) const;

    void
    DL_invM (SnVector &dst,
             const SnVector &src) const;

    /**
     * @brief parallel version of DL_invMS
     */
    void
    ParDL_invMS (BlockVector<double> &dst,
                 const BlockVector<double> &src) const;

    /**
     * @brief parallel version of DL_invMS
     */
    void
    ParDL_invM (BlockVector<double> &dst,
                const BlockVector<double> &src) const;

    /** @brief reference to the State object provided as argument. */
    Forest::State<dim> & mp_state;
    /** @brief reference to solver used internally to provide the .solve() method . */
    mutable EqSolverWG<dim> m_solver;
    /** @brief Do we want to do something else than matrix free? */
    bool m_matrix_free;
    /** @brief Constant reference to the quadrature object. */
    std::shared_ptr<QuadratureBase<dim> > m_ord;
    /** @brief Transport operator. */
    std::shared_ptr<TransOp<dim> > m_L0;
    /** @brief Transport bcs. */
    std::shared_ptr<BoundaryValues<dim> > m_BV;
    /** @brief Fission spectrum. */
    std::shared_ptr<FissionSpectrum<dim> > m_Chi;
    /** @brief Fission source. */
    std::shared_ptr<FissionSource<dim> > m_XF;
    /** @brief Scattering source.  */
    std::shared_ptr<ScatteringSource<dim> > m_S;
    /** @brief Moments-to-Directions operator. */
    mutable MomentsToDirections<dim> m_M;
    /** @brief Directions-to-Moments operator. */
    mutable DirectionsToMoments<dim> m_D;

    /** @brief Preconditioner. */
    mutable std::shared_ptr<Equation<dim> > m_prec;

    mutable unsigned int m_g;
    mutable bool m_external_preconditioner;
    /**
     * @brief Vector containing the Spherical Harmonics moments of the flux.
     * @details The size of this vector is \f$ I \times L \f$, i.e., the number
     * of elements for the spatial mesh times the number of Legendre moments.
     */
    mutable BlockVector<double> m_v0, m_v1;
    mutable Vector<double> m_v2, m_v3;
    mutable BlockVector<double> m_pv2, m_pv3;

    mutable std::vector<unsigned int> m_angular_sizes;
    mutable std::vector<unsigned int> m_flux_bv_sizes;
    mutable std::vector<unsigned int> m_fission_density_sizes;

  public:
    virtual bool no_upscatt(const unsigned int g) const;
  };

} // end of namespace Forest

#endif /* FOREST_EQ_TRANSPORT_H */
