/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eigen_prob.h
 * @brief   EigenProb class template declarations
 */

#ifndef FOREST_EIGEN_PROB_H
#define FOREST_EIGEN_PROB_H

#include "algebra/sn_vector_fwd.h"
#include "algebra/equation.h"
#include "algebra/solver_eq.h"
#include "neutronics/state.h"

#include <string>
#include <vector>
#include <memory>

#ifndef DOXYGEN
namespace dealii { template <typename > class BlockVector; }
namespace dealii { template <typename > class Vector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
    @class EigenProb
    @ingroup ForestAlgebra
    @brief Eigenvalue problem; base class.

    @details Base class to construct the derived eigenproblems. This class
    specifies the contract to be satisfied by the derived class. In
    general we are solving a generalized eigenvalue problem, but we want
    to rearrange the different parts to end up with a system of the form

    \f{align}{
     A x  = \lambda x ,
    \f}

    ----------

    For instance, we can think about the generalized eigenvalue problem
    derived from the neutron transport equation

    \f{align}{
     (DLM + S) \phi  = \frac{1}{\lambda} X F^{T} \phi ,
    \f}

    This problem can be rearrange in different ways as a standard
    eigenvalue problem:

      - We can invert the right hand side (it must always have an inverse),
        so we then solve for the problem

        \f{align}{
        A \phi = \lambda \phi ,
        \quad \text{where}
        \quad A :=  (DLM + S)^{-1} X F^{T} .
        \f}


      - The problem can also be rearranged as

        \f{align}{
        A f = \lambda f ,
        \quad \text{where}
        \quad f := F^{T} \phi
        \quad \text{and}
        \quad A := F^{T} (DLM + S)^{-1} X .
        \f}

    Here we also store some temporary vectors that we need when performing
    the different operations. This temporary vectors would be for
      - @p m_phi_tmp the scalar flux + boundary values
      - @p m_fission_density_tmp
      - @p m_fission_source_tmp

    @note In the particular case when the operator M is built by the product
    of two reduced rank matrices, as for example when we multiply the fission
    spectrum times the fission cross section, \f$ M = X F^{T} \f$, the
    eigenvalue can be recast in terms of the fission density
    \f$ f = F^{T}\phi \f$, to obtain the eigenproblem where the matrix
    \f$ A  := F^{T} (DLM + S)^{-1} X \f$ would have the
    dimensions of one energy group. This can potentially increase the
    efficiency of the Krylov solver by reducing the dimensionality of the
    search space for the solution.

    Reference for Jacobian free:

    > A comparison of acceleration methods for solving the neutron transport
    > k-eigenvalue problem;
    > Jeffrey Willert , H. Park, D.A. Knoll;
    > Journal of Computational Physics 274 (2014) 681–694
   */
  template <int dim>
  class EigenProb
  {
  public:


    /**
     * @brief Constructor of the eigenproblem
     * @param state
     * @param equation
     */
    EigenProb (State<dim> & state,
               std::shared_ptr<Equation<dim> > & equation)
        : mp_state (state),
          mp_equation (equation),
          m_solver (mp_state, mp_equation)
    {
      const unsigned int n_groups = mp_state.get_n_groups ();
      m_phi_tmp.reinit (n_groups, mp_equation->get_flux_bv_sizes ());
      m_fission_source_tmp.reinit (n_groups,
          mp_equation->get_fission_density_sizes ());
      m_fission_density_tmp.reinit (1, mp_equation->get_fission_density_sizes ());
    }

    /** @brief Destructor. */
    virtual
    ~EigenProb ()
    {
    }

    /** @brief Pure virtual function defining the matrix vector product. */
    virtual void
    vmult (SnVector &dst,
           const SnVector &src) const = 0;

    /** @brief The transpose operator. */
    virtual void
    Tvmult (SnVector & dst,
            const SnVector & src) const = 0;

    /** @brief Resize the vector with the size of the solution. */
    virtual void
    set_solution_sizes (SnVector & vec) const = 0;

    /** @brief Set the initial guess. */
    virtual void
    set_initial_guess () const = 0;

    /** @brief Access function for the initial guess */
    virtual SnVector &
    get_initial_guess () const = 0;

    /**
     @brief We update the value of the vectors using the solution.

     @details If we have used an Arnoldi method for the eigenvalue problem,
     the intermediary vectors would be vectors of the krylov subspace.
     In order for the intermediary vectors to have a physical meaning
     we must perform one extra iteration using the solution vector. It have been
     implemented only in derived classes.
    */
    virtual void
    update_vectors(const SnVector & vec, const double keff) const = 0;

    /** @brief Access function for the fission density */
    SnVector
    get_fission_density () const
    {
      return m_fission_density_tmp;
    }

    /** @brief Access function for phi */
    SnVector
    get_phi () const
    {
      return m_phi_tmp;
    }

    /**
      @brief Access function for psi.
      @details We do not have the angular neutron flux stored
      because of the large amount of memory resources that it might
      require. If anyway we want to access it, we will calculate
      it using the scalar neutron flux and the "incoming angular flux"
      at the boundaries, and the container vector must be provided
      from outside even if it will be resized inside this function.
      This operation requires one full transport sweep.
      */
    void
    get_psi (SnVector & psi, const double keff) const
    {
      const unsigned int n_groups = mp_state.get_n_groups ();
      psi.reinit (n_groups, mp_equation->get_angular_sizes ());
      mp_equation->get_angular (psi, m_phi_tmp, keff);
    }

    /** @brief Access to m_solver.get_sweeps_per_group() */
    double
    get_sweeps_per_group () const
    {
      return m_solver.get_sweeps_per_group ();
    };

    /** @brief Access to m_solver.get_sweeps_per_group() */
    unsigned int
    get_total_sweeps () const
    {
      return m_solver.get_total_sweeps();
    }

    /** @brief Memory consumption of object */
    virtual
    double
    memory_consumption () const = 0;

  protected:

    /** @brief reference to the State object provided as argument. */
    Forest::State<dim> & mp_state;
    /** @brief Constant reference to the equation object. */
    std::shared_ptr<Equation<dim> > & mp_equation;
    /** @brief Container for the solver. */
    SolverEq<dim> m_solver;

    /** @brief Container for the flux + bv. */
    mutable SnVector m_phi_tmp;
    /** @brief Container for the fission density. */
    mutable SnVector m_fission_density_tmp;
    /** @brief Container for the fission source. */
    mutable SnVector m_fission_source_tmp;
  };

} // end of namespace Forest

#endif /* FOREST_EIGEN_PROB_H */
