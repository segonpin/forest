/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eigen_solver.h
 * @brief   EigenSolver class template declarations
 */

#ifndef FOREST_EIGEN_SOLVER_H
#define FOREST_EIGEN_SOLVER_H

#include "neutronics/state_fwd.h"       // For class State&
#include "algebra/eigen_prob.h"
#include "neutronics/extract_bv.h"
#include "angle/moments_to_directions.h"
#include "angle/directions_to_moments.h"
#include "algebra/equation.h"
#include "algebra/eigen_prob_factory.h"

#include <memory>

#ifndef DOXYGEN
namespace Forest { template <int dim> class ProbGeom; }
namespace Forest { class Input; }
#endif

#include "algebra/sn_vector_fwd.h"
#include "algebra/equation_fwd.h"       // For class Equation&


namespace Forest
{

  /**  
    @class EigenSolver
    @ingroup ForestAlgebra
    @brief Eigenvalue problem solver

    @details An object @p eigen_prob of class EigenProb
    is provided that defines a standard eigenvalue problem.
    This eigenvalue problem is then solved with any of the methods
    in this class.
   */
  template <int dim>
  class EigenSolver
  {
  public:

    /**
     @brief Use the state object to initialize basic information.
     */
    EigenSolver (State<dim> &state,
                 std::shared_ptr<EigenProb<dim> > & eigen_prob);

  public:

    /**
     @brief This method runs the power iteration function.
     */
    void
    solve_eigenproblem ();

    /**
      @brief Power iteration method for the eigenvalue problem.

      @details Given the problem

      \f{align*}{
      A x = \lambda x
      \f}

      we use the power iteration to approximate the pair \f$ (\lambda, x) \f$.

        - We start with an initial guess \f$ (\lambda_0, x_0) \f$,
          and the first step is to calculate the matrix vector product
          \f{align*}{
          & y = A x_i
          \f}
          and update the eigen-pair
          \f{align*}{
          & \lambda_{i+1} = y^{t} \cdot x_i \ , \quad
          x_{i+1} = \frac{y}{\lambda_{i+1}}
          \f}

        - Now we estimate the different errors
          \f{align*}{
          & e_{\lambda} = \frac{\| \lambda_{i+1} - \lambda_{i} \| }{\lambda_{i}} \ , \quad
          e_{x} = \| x_{i+1} - x_{i} \|
          \f}
          and decide the stopping criteria. If this criteria is not fulfilled,
          we increase the iteration index and start again.
     */
    void
    power_iteration (SnVector &phi, double &keff);


    /**
     * @brief Solve the eigenvalue problem with SLEPc, by default using the Krylov-Schur
     * method.
     *
     * @details SLEPc need its own PETSc vector so a full vector copy must be performed.
     *  This copy was tested to cost less than 0.1% of the total computational time.
     */
    void
    slepc (const unsigned int n_eigenvalues,
           const unsigned int n_cv,
           SnVector &phi,
           double &keff);


  private:

    /**
     * @brief Message about convergence.
     * @param keff_converged
     */
    void
    convergence_msg (const bool keff_converged);

    /** @brief returns the memory consumption of the objects in this class. */
    double
    memory_consumption () const;

    /**
     * @brief reference to the State object provided as argument.
     * @details It is not constant because we have to change the solution, which
     * is inside the state object.
     */
    Forest::State<dim> & mp_state;

    /** @brief Container for the eigenproblem. */
    std::shared_ptr<EigenProb<dim> > m_eigen_prob;

    /** @brief Container for the eigenvalue. */
    double m_keff;

    /** @brief Container for the fission source. */
    SnVector m_tmp_old, m_tmp_new;

    unsigned int m_max_it;
    double m_tol;


  };

} // end of namespace Forest
#endif /* FOREST_EIGEN_SOLVER_H */
