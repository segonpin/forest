/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eigen_solver.cc
 * @brief  Implementation of class template EigenSolver
 */

#include "neutronics/state.h"
#include "algebra/sn_vector.h"
#include "input/input.h"
#include "algebra/solver_eq.h"
#include "utils/forest_utils_logstream.h"
#include "utils/forest_utils_memory.h"
#include "utils/forest_utils_timer.h"
#include "utils/forest_utils_base.h"
#include "algebra/eigen_solver.h"
#include "algebra/eigen_prob_factory.h"
#include "algebra/eigen_solver_slepc.h"

#include "deal.II/base/timer.h"
#include "deal.II/lac/slepc_solver.h"

#include "deal.II/base/mpi.h"

#include <memory>
#include <cmath>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  EigenSolver<dim>::EigenSolver (State<dim> & state,
                                 std::shared_ptr<EigenProb<dim> > & eigen_prob)
      : mp_state (state),
        m_eigen_prob (eigen_prob),
        m_keff (1.),
        m_max_it (mp_state.mp_data.mp_settings.get_eig_max_it ()),
        m_tol (mp_state.mp_data.mp_settings.get_eig_tol ())
  {
    bool verbose = false;
    if (verbose)
    {
      log_print_text ("Setup and assemble systems in EigenSolver", 1);
    }

    // Set the size of the auxiliary vectors.
    m_eigen_prob->set_solution_sizes (m_tmp_new);
    m_eigen_prob->set_solution_sizes (m_tmp_old);

    // Printing Stats. for the Memory used by the auxiliary vectors.
    extern StatsMemory memory;
    StatsMemory::instance ().add ("EigenSolver", memory_consumption ());
    StatsMemory::instance ().print_detailed ();
  }

  template <int dim>
  double
  EigenSolver<dim>::memory_consumption () const
  {
    double memory = 0;
    memory += m_tmp_new.memory_consumption () + m_tmp_old.memory_consumption ();
    return memory;
  }

  template <int dim>
  void
  EigenSolver<dim>::power_iteration (SnVector & vec,
                                     double & keff)
  {
    Timer timer;
    timer.restart ();

    // Declare variables for the keff.
    double keff_old, err_keff, err_vec;

    // Get the initial guess.
    m_tmp_new = vec;
    m_tmp_new *= 1.0 / m_tmp_new.l2_norm ();

    log_print_line ();
    bool keff_converged = false;
    for (unsigned int iter = 0; iter < m_max_it; ++iter)
    {
      // Move new (or guess) to old solution.
      m_tmp_old = m_tmp_new;
      keff_old = keff;

      // Apply the operator y = A x.
      m_eigen_prob->vmult (m_tmp_new, m_tmp_old);

      // Update the k-effective and the eigenvector.
      keff = m_tmp_new * m_tmp_old; //this->update_keff (m_tmp_new, m_tmp_old);
      m_tmp_new *= 1.0 / m_tmp_new.l2_norm ();
      m_tmp_old.add (-1.0, m_tmp_new);

      // estimate the errors
      err_keff = std::abs (keff - keff_old) / std::abs (keff_old);
      err_vec = m_tmp_old.l2_norm ();
      const double error_comp = err_vec + err_keff;

      // Print information about the current mp_state of the convergence.
      const double spg = m_eigen_prob->get_sweeps_per_group ();
      keff_converged = log_table_keff_wgit (iter, keff, spg, error_comp, m_tol);
      // keff_converged = log_table_keff (iter, keff, err, m_tol);

      // Stop the loop if convergence has been achieved.
      if (keff_converged)
      {
        log_print_text (
            "Total number of iterations: " + std::to_string (iter + 1), 1);
        break;
      }
    }

    // Information about the number of sweeps and the convergence achieved.
    // Move the solution to the container to be returned.
    const unsigned int n_total_sweeps = m_eigen_prob->get_total_sweeps ();
    log_print_text (
        "Total number of sweeps: " + std::to_string (n_total_sweeps), 1);
    this->convergence_msg (keff_converged);
    vec = m_tmp_new;

    timer.stop ();
    StatsTimer::instance ().add ("Power iteration", timer.wall_time ());
  }

  /**
   * Call SLEPc eigenvalue problem solve to solve.
   * TODO Change interaction to std::vector<SnVector> and std::vector<double>
   */
  template <int dim>
  void
  EigenSolver<dim>::slepc (const unsigned int n_eigenvalues,
                           const unsigned int n_cv,
                           SnVector & vec,
                           double & keff)
  {
    Timer timer;
    timer.restart ();



    // Initializes vectors.
    std::vector<double> eigenvalues (n_eigenvalues);
    std::vector<SnVector> psi (n_eigenvalues, vec);
    psi[0] = vec;

    EigenSolverSLEPC<dim> slepc_solver (m_eigen_prob, n_eigenvalues, n_cv);
    slepc_solver.tol_eps = m_tol;
    int error_code = slepc_solver.solve (eigenvalues, psi);
    bool keff_converged = (error_code == 0);

    // Information regarding the number of iterations and sweeps.
    log_print_line ();
    log_print_text (
        "Total number of iterations: " + std::to_string (
            slepc_solver.n_matvec_multiplications), 1);
    const unsigned int n_total_sweeps = m_eigen_prob->get_total_sweeps ();

    log_print_text ("Total number of sweeps: " + num_to_string (n_total_sweeps),1);
    log_print_text ("Krylov subspace dimension : " + num_to_string (n_cv), 1);
    log_print_text ("Number of eigenvalues: " + num_to_string (n_eigenvalues), 1);

    // Print Results in log.
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      log_print_text (
          "k" + num_to_string (eig + 1) + " = "
          + num_to_string (eigenvalues[eig]), 1);

    // Copy back vectors
    keff = eigenvalues[0];
    vec = psi[0];

    // Information about the convergence achieved by the loop.
    this->convergence_msg (keff_converged);

    // Stop timer and report it.
    timer.stop ();
    StatsTimer::instance ().add ("Krylov iteration", timer.wall_time ());
  }

  /*template <int dim>
   void
   EigenSolver<dim>::jacobi_davidson (SnVector & phi,
   double & keff)
   {
   // Start timer
   Timer timer;
   timer.restart ();

   // Print line
   log_print_line ();
   log_print_text ("Solving with Jacobi-Davidson", 1);

   // Defining things
   const PETScWrappers::MatrixFree A;
   std::vector< PetscScalar > lambda;
   std::vector< PetscScalar > eigenvectors;
   const unsigned int n_eigenpairs=1;

   // Solve
   SolverControl solver_control (1000, 1e-9);
   SLEPcWrappers::SolverJacobiDavidson system (solver_control);
   system.solve (A, lambda, eigenvectors, n_eigenpairs);

   // Recover
   phi = eigenvectors[0];
   keff = lambda[0];

   // Print line
   log_print_line ();

   // Stop timer and report it.
   timer.stop ();
   StatsTimer::instance ().add ("Jacobi-Davidson", timer.wall_time ());
   }
   */

  template <int dim>
  void
  EigenSolver<dim>::convergence_msg (const bool keff_converged)
  {
    if (keff_converged)
      log_print_text ("Keff converged for given tolerance.", 1);
    else
      log_print_text (
          "Keff NOT converged. Increase 'eig_max_it' or reduce 'eig_tol'.", 1);
  }

  template <int dim>
  void
  EigenSolver<dim>::solve_eigenproblem ()
  {
    // Initialize the keff to 1. (if no better guess).
    m_keff = 1.;
    // Initialize fission density with flux 1.0 (if no better guess).
    m_eigen_prob->set_initial_guess ();
    SnVector & solution = m_eigen_prob->get_initial_guess ();
    //m_keff = solution.l2_norm();
    //solution *=1./solution.l2_norm();

    const unsigned int n_eigenvalues =
        mp_state.mp_data.mp_settings.get_n_eigenvalues ();
    const unsigned int n_cv =
        mp_state.mp_data.mp_settings.get_n_column_vectors();


    // Call the eigenvalue solver.
    if (mp_state.mp_data.mp_settings.use_power_iteration ())
      this->power_iteration (solution, m_keff);
    else
      this->slepc (n_eigenvalues, n_cv, solution, m_keff);
    m_eigen_prob->update_vectors (solution, m_keff);
    log_print_text ("After updating the vectors", 1);

    // Set vector for scalar flux together with the boundary values.
    SnVector & phi = mp_state.m_scalar_flux;
    phi = m_eigen_prob->get_phi ();

    // Store the solution back in the state class.
    mp_state.set_keff (m_keff);

    // TODO
    if (mp_state.mp_data.mp_settings.use_power_iteration ())
    {
      // Apply one sweep to recover the angular flux.
      SnVector & psi = mp_state.m_psi;
      m_eigen_prob->get_psi (psi, m_keff);
    }
    else
      log_print_text ("If slepc we jump to avoid crashing", 1);
  }

  template class EigenSolver<1> ;
  template class EigenSolver<2> ;
  template class EigenSolver<3> ;

} // end of namespace Forest
