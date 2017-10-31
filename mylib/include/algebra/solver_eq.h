/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/solver_eq.h
 * @brief  SolverEq, EqSolverWG and EqPrec class template declarations
 */

#ifndef FOREST_SOLVER_EQ_H
#define FOREST_SOLVER_EQ_H

#include "algebra/sn_vector.h"          // for SnVector
#include "algebra/equation_fwd.h"       // for Equation &
#include "neutronics/state_fwd.h"       // For class State&

#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/lac/vector.h"         // for Vector
#include "deal.II/lac/vector_memory.h"  // for GrowingVectorMemory

#include <vector>
#include <memory>

namespace Forest
{
  using namespace dealii;

  /**
   * @class EqPrec
   * @brief Preconditioner structure to wrap the .prec() member of equation inside with a vmult use by deal.II
   * @ingroup ForestAlgebra
   */
  template <int dim>
  class EqPrec
  {
  public:
    EqPrec (Equation<dim> & equation)
        : m_equation (equation)
    {
    }
    void
    vmult (BlockVector<double> & phi,
           const BlockVector<double> & source) const
    {
      m_equation.prec (phi, source);
    }
  private:
    const Equation<dim> & m_equation;
  };


  template <int dim>
  class Multigroup
  {
  public:
    Multigroup (std::shared_ptr<Equation<dim> > & equation)
        : m_equation (equation)
    {
    }

    void
    vmult (SnVector & dst,
           const SnVector & src) const
    {
      dst = 0;
      m_equation->vmult (dst, src);
    }

    void
    prec(SnVector & phi_new,
         const SnVector & mg_rhs) const
    {
      //prec_gs(phi_new, mg_rhs);
      prec_jac(phi_new, mg_rhs);
      //phi_new = mg_rhs;
    }

    void
    prec_gs(SnVector & phi_new,
        const SnVector & mg_rhs) const
    {
      BlockVector<double> m_outer_source;
      BlockVector<double> m_rhs, m_sol;

      // Loop over the energy groups to solve the source iteration per group.
      for (unsigned int g = 0; g < m_equation->get_n_groups (); ++g)
      {
        // Ther within group source starts with the fission using the old flux
        m_outer_source = mg_rhs.m_data[g];
        m_outer_source = 0; // this is just to make the zero with the right size
        m_rhs = m_outer_source; // also to make it zero

        // Add down-scatt. with new flux (if Gauss-Seidel)
        m_equation->down_add (m_outer_source, phi_new, g);
        // We tell the energy group to the matrix vector product.
        m_equation->set_group (g);
        // First we prepare the right hand side
        m_equation->set_rhs (m_rhs, m_outer_source);
        // we add one sweep that is needed to generate the rhs
        m_rhs.sadd (-1., mg_rhs.m_data[g]);

        // Put the flux and the boundary conditions in one vector.
        m_sol = phi_new.m_data[g];

        // Use Krylov method or Richardson iteration.
        // const double new_tol = m_tol_i / outer_error;
        // const unsigned int new_max_it = m_max_it_i;
        /// @todo add new_tol, new_max_it?
        deallog << "g = " << g << "; ";
        m_equation->solve (m_sol, m_rhs);

        // Extract the solution and the boundary conditions from the one vector.
        phi_new.m_data[g] = m_sol;

        // Accumulate the within group iterations to the total iterations.
        // m_total_sweeps_this_iteration += m_equation->get_tot_wg_it ();

      }
    }

    void
    prec_jac(SnVector & phi_new,
        const SnVector & mg_rhs) const
    {
      BlockVector<double> m_rhs, m_sol;

      // Loop over the energy groups to solve the source iteration per group.
      for (unsigned int g = 0; g < m_equation->get_n_groups (); ++g)
      {
        m_equation->set_group (g);
        deallog << "g = " << g << "; ";
        // we add one sweep that is needed to generate the rhs
        m_rhs = mg_rhs.m_data[g];
        // Put the flux and the boundary conditions in one vector.
        m_sol = phi_new.m_data[g];

        // Use Krylov method or Richardson iteration.
        // const double new_tol = m_tol_i / outer_error;
        // const unsigned int new_max_it = m_max_it_i;
        /// @todo add new_tol, new_max_it?
        m_equation->solve (m_sol, m_rhs);

        // Extract the solution and the boundary conditions from the one vector.
        phi_new.m_data[g] = m_sol;

        // Accumulate the within group iterations to the total iterations.
        // m_total_sweeps_this_iteration += m_equation->get_tot_wg_it ();
      }
    }
  private:
    std::shared_ptr<Equation<dim> > & m_equation;
    //const Equation<dim> & m_equation;
  };


  template <int dim>
  class MGPrec
  {
  public:
    MGPrec (Multigroup<dim> & multigroup_equation)
        : m_equation (multigroup_equation)
    {
    }
    void
    vmult (SnVector & phi,
           const SnVector & source) const
    {
      m_equation.prec (phi, source);
    }
  private:
    const Multigroup<dim> & m_equation;
  };

  /**
   * @brief Class implementing the Solver.
   * @ingroup ForestAlgebra
   * @details This class implements the solver for the source problem
   * working at two different levels, the within-group source problem
   * and the multi-group problem, depending on the type of the vector provided.
   *
   * We have to solve a source problem in a multi-group level,
   * and this is approximated by working with single-group solvers.
   *
   * # Multi-Group (Gauss-Seidel)
   *
   * The original problem looks as follows (integro-differential operator)
   *
   * \f{align*}{
   * (L+MSD)\Psi = Q
   * \f}
   *
   * This system has a block structure defined by the energy discretization,
   * where the streaming operator \f$ L \f$ is block diagonal (no coupling
   * between different groups), and operator \f$ S \f$ is "almost" lower
   * triangular (up-scattering is usually smaller than down-scattering). A
   * reasonable approach would be then to solve the systems as if it would be
   * lower triangular, and "lagger" the up-scattering to be multiplied to a
   * guess. To do so, we start solving the system for the first energy group,
   * and then applying the Gauss-Seidel method as follows
   *
   * \f{align*}{
   * (L+MS_{LD}D)\Psi^{n+1} = M S_{U} D \Psi^{n} + Q
   * \f}
   *
   * But, in order to improve its solvability the system is rewritten as
   *
   * \f{align*}{
   * (I+L^{-1}MS_{LD}D)\Psi^{n+1} =  L^{-1}( M S_{U} D \Psi^{n} + Q )
   * \f}
   *
   * Or we can simplify the notation by writing
   *
   * \f{align*}{
   * A x = b
   * \f}
   *
   * where the matrix \f$ A \f$ (approximating an integral operator)
   * and the source \f$ b \f$ are defined as
   *
   * \f{align*}{
   * A := (I-L^{-1}MS_{LD}D) ,
   * \quad \text{and} \quad
   * b := L^{-1}\hat{Q} := L^{-1} (M S_{U} D \Psi^{n} + Q)
   * \f}
   *
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
   * ## Neglecting Reflecting Boundary Conditions.
   *
   * We use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality, and then update the solution inverting the
   * streaming term (what should be easier)
   *
   * \f{align*}{
   * & \Psi^{l+1} = L^{-1}(MS \Phi^{l} + Q) \\
   * & \Phi^{l} = D \Psi^{l}
   * \f}
   *
   * This expression, ignoring boundary conditions, can be simplified further
   * by combining the two equations as follows
   *
   * \f{align*}{
   * & (I - D L^{-1} MS) \Phi = D b
   * \f}
   *
   *
   * ## Including Reflecting Boundary Conditions
   *
   * Here we assume "lagged" the scattering term to a previous iteration step,
   * as before, but also "lagged" the effect of the boundary conditions, and
   * we use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality. Thus we obtain
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
   * D L_{0}^{-1} Q \\
   * P^{T}L_{0}^{-1} Q
   * \end{array} \right]
   * \f}
   *
   * This can be rewritten as
   *
   * \f{align*}{
   * & (\hat{I} - \hat{D} \hat{L}^{-1} \hat{M}\hat{S}) \hat{\Phi} = \hat{b}
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
   * @todo Hide all the discussion related to the particular problem from this
   * class and move it to the Equation class. Here we are solving a linear
   * system \f$ A x = b \f$ .
   */
  template <int dim>
  class EqSolverWG
  {
  public:

    /** @todo Document me */
    EqSolverWG (State<dim> & state,
                Equation<dim> & equation);

    /** @todo Document me */
    void
    solve (BlockVector<double>  & sol,
           const BlockVector<double> & rhs,
           const double tol_default = -1.,
           const unsigned int max_it_default = 0) const;

    /**
      @brief Richardson method for solving a linear system
      @details
      phi_new = [D*L_inv*M*S]*phi_old
              = [I-(I-D*L_inv*M*S)]*phi_old
              = [I-Eq.vmult]*phi_old ;
     */
    void
    richardson (BlockVector<double> & sol,
                const BlockVector<double> & rhs,
                const double tol,
                const unsigned int max_it) const;
  public:

    /** @brief Return the total number of sweeps within energy groups. */
    unsigned int
    get_tot_wg_it () const
    {
      return m_total_wg_it;
    }

    /** @brief Return the memory consumption of the vectors of this class. */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      memory_consumption += m_vector_memory.memory_consumption();
      return memory_consumption;
    }

  private:

    /** @brief Reference pointer to the state vector. */
    State<dim> & mp_state;
    Equation<dim> & mp_equation;
    EqPrec<dim> m_preconditioner;

    /** @brief Memory for GMRES (within-group). */
    mutable GrowingVectorMemory<BlockVector<double> > m_vector_memory;

    /** @brief Maximum iterations of the solvers (within-group). */
    const unsigned int m_max_it_i;
    /** @brief Tolerance criteria for the solvers (within-group). */
    const double m_tol_i;
    /** @brief Flag to use the Krylov solver (within-group). */
    bool m_use_krylov;
    /** @brief Flag to use the Richardson iteration (within-group). */
    bool m_use_richardson;
    /** @brief Counting the sweeps (within-group). */
    mutable unsigned int m_total_wg_it;
  };

  /**
   * @brief Class implementing the Solver.
   * @ingroup ForestAlgebra
   * @details This class implements the solver for the source problem
   * working at two different levels, the within-group source problem
   * and the multi-group problem, depending on the type of the vector provided.
   *
   * We have to solve a source problem in a multi-group level,
   * and this is approximated by working with single-group solvers.
   *
   * # Multi-Group (Gauss-Seidel)
   *
   * The original problem looks as follows (integro-differential operator)
   *
   * \f{align*}{
   * (L+MSD)\Psi = Q
   * \f}
   *
   * This system has a block structure defined by the energy discretization,
   * where the streaming operator \f$ L \f$ is block diagonal (no coupling
   * between different groups), and operator \f$ S \f$ is "almost" lower
   * triangular (up-scattering is usually smaller than down-scattering). A
   * reasonable approach would be then to solve the systems as if it would be
   * lower triangular, and "lagger" the up-scattering to be multiplied to a
   * guess. To do so, we start solving the system for the first energy group,
   * and then applying the Gauss-Seidel method as follows
   *
   * \f{align*}{
   * (L+MS_{LD}D)\Psi^{n+1} = M S_{U} D \Psi^{n} + Q
   * \f}
   *
   * But, in order to improve its solvability the system is rewritten as
   *
   * \f{align*}{
   * (I+L^{-1}MS_{LD}D)\Psi^{n+1} =  L^{-1}( M S_{U} D \Psi^{n} + Q )
   * \f}
   *
   * Or we can simplify the notation by writing
   *
   * \f{align*}{
   * A x = b
   * \f}
   *
   * where the matrix \f$ A \f$ (approximating an integral operator)
   * and the source \f$ b \f$ are defined as
   *
   * \f{align*}{
   * A := (I-L^{-1}MS_{LD}D) ,
   * \quad \text{and} \quad
   * b := L^{-1}\hat{Q} := L^{-1} (M S_{U} D \Psi^{n} + Q)
   * \f}
   *
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
   * ## Neglecting Reflecting Boundary Conditions.
   *
   * We use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality, and then update the solution inverting the
   * streaming term (what should be easier)
   *
   * \f{align*}{
   * & \Psi^{l+1} = L^{-1}(MS \Phi^{l} + Q) \\
   * & \Phi^{l} = D \Psi^{l}
   * \f}
   *
   * This expression, ignoring boundary conditions, can be simplified further
   * by combining the two equations as follows
   *
   * \f{align*}{
   * & (I - D L^{-1} MS) \Phi = D b
   * \f}
   *
   *
   * ## Including Reflecting Boundary Conditions
   *
   * Here we assume "lagged" the scattering term to a previous iteration step,
   * as before, but also "lagged" the effect of the boundary conditions, and
   * we use the flux in the previous iteration collapsed with the operator D
   * to reduce its dimensionality. Thus we obtain
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
   * D L_{0}^{-1} Q \\
   * P^{T}L_{0}^{-1} Q
   * \end{array} \right]
   * \f}
   *
   * This can be rewritten as
   *
   * \f{align*}{
   * & (\hat{I} - \hat{D} \hat{L}^{-1} \hat{M}\hat{S}) \hat{\Phi} = \hat{b}
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
   * @todo Hide all the discussion related to the particular problem from this
   * class and move it to the Equation class. Here we are solving a linear
   * system \f$ A x = b \f$ .
   */
  template <int dim>
  class SolverEq
  {
  public:

    /** @todo Document me */
    SolverEq (State<dim> & state,
              std::shared_ptr<Equation<dim> > & equation);

    /** @todo Document me */
    void
    iterate (std::shared_ptr<Equation<dim> > & equation,
             SnVector & phi_new,
             const SnVector & fis_source_old,
             const double outer_error = 1.) const;

    /** @todo Document me */
    void
    gauss_seidel_LDU (std::shared_ptr<Equation<dim> > & equation,
             SnVector & phi_new,
             const SnVector & phi_old,
             const SnVector & fis_source_old) const;

  public:

    /** @brief Return the average number of sweeps per energy groups. */
    double
    get_sweeps_per_group () const;

    /** @brief Return the total number of sweeps for all energy groups. */
    unsigned int
    get_total_sweeps () const
    {
      return m_total_sweeps;
    }

    /** @brief Return the memory consumption of the vectors of this class. */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      memory_consumption +=
          + m_rhs.memory_consumption ()
          + m_sn_phi.memory_consumption ()
          + m_tmp_source.memory_consumption ()
          + m_aux_vec.memory_consumption ()
          + m_outer_source.memory_consumption ()
          + m_sol.memory_consumption ()
          + m_phi_old.memory_consumption ()
          + m_sol_old.memory_consumption ();
      return memory_consumption;
    }

  private:

    /** @brief Reference pointer to the state vector. */
    State<dim> & mp_state;

    /** @brief Memory for GMRES (within-group). */
    mutable GrowingVectorMemory<BlockVector<double> > m_vector_memory;
    mutable GrowingVectorMemory<SnVector > m_mg_vector_memory;

    /** @brief Maximum iterations of the solvers (within-group). */
    const unsigned int m_max_it_i;
    /** @brief Tolerance criteria for the solvers (within-group). */
    const double m_tol_i;

    /** @brief Flag to use the Krylov solver (within-group). */
    bool m_use_krylov;
    /** @brief Flag to use the Gauss-Seidel scheme (full-energy). */
    bool m_use_gauss_seidel;


    /** These constants are used to initialize some vectors. */
    unsigned int m_n_leg_mom;
    unsigned int m_n_elem;
    unsigned int m_n_groups;

    /** @brief Accumulating the sweeps for all groups (full-energy). */
    mutable unsigned int m_total_sweeps_this_iteration;
    mutable unsigned int m_total_sweeps;
    /** @brief Auxiliary vectors. */
    mutable Vector<double> m_sn_phi, m_tmp_source;
    /** @brief Auxiliary vectors. */
    mutable BlockVector<double> m_aux_vec;

    /** @brief Size of the boundary conditions per direction. */
    mutable std::vector<unsigned int> m_flux_bv_sizes;

    /** @brief Declaring the rhs vector. */
    mutable BlockVector<double> m_rhs;
    mutable SnVector m_mgrhs;
    /** @brief Declaring the solution vector. */
    mutable BlockVector<double> m_sol;
    /** @brief Container for the old solution. */
    mutable BlockVector<double> m_sol_old;

    mutable BlockVector<double> m_outer_source;

    /** @brief Container for the scalar flux. */
    mutable SnVector m_phi_old;
  };

} // end of namespace Forest

#endif /* FOREST_SOLVER_EQ_H */
