/**
 * @author  Antoni Vidal. UPV, 2017.
 * @file    algebra/eigen_solver_slepc.cc
 * @date    Modified on 23/01/2017
 * @brief   EigenSolver class to solve the function with SLEPCS
 */

#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/lapack_support.h>

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

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>

#include <stdlib.h>
#include <stdio.h>

#include "algebra/eigen_solver_slepc.h"

#include <petscviewer.h>

namespace Forest
{

  /**
   *
   */
  template <int dim>
  EigenSolverSLEPC<dim>::EigenSolverSLEPC (std::shared_ptr<EigenProb<dim> > _eigen_prob,
                                           const unsigned int _n_eigenvalues,
                                           const unsigned int _n_cv)
      : m_eigen_prob (_eigen_prob),
        n_eigenvalues (_n_eigenvalues),
        n_cv (_n_cv)
  {
    n_iterations = 0;
    n_matvec_multiplications = 0;
    n_groups = 0;

    // Default values
    tol_eps = 1e-7;
    max_iterations_eps = 200;

    // Silly initializations
    shell_mat = NULL;
    eps = NULL;

    tmp_vector = NULL;

    // -ksp_monitor
    // -ksp_converged_reason

    return;
  }

  /**
   * Destroy The ksp, eps and shell objects.
   */
  template <int dim>
  EigenSolverSLEPC<dim>::~EigenSolverSLEPC ()
  {
    EPSDestroy (&eps);
    MatDestroy (&shell_mat);
    return;
  }

  /**
   *
   */
  template <int dim>
  int
  EigenSolverSLEPC<dim>::solve (std::vector<double>& eigenvalues,
                                std::vector<SnVector>& psi)
  {
    PetscErrorCode ierr;

    // Initialize EPS
    ierr = EPSCreate (PETSC_COMM_WORLD, &eps);
    // TODO Set ncv as as and input parameter
    // Default 2*n_ev + 1
    ierr = EPSSetDimensions (eps, n_eigenvalues, n_cv, PETSC_DEFAULT);

    // Copy initial vector
    sn_vector_tmp_src = psi[0];
    n_groups = psi[0].size0 ();

    // Calculate n_dofs
    int n_dofs = 0;
    for (unsigned int g = 0; g < n_groups; ++g)
      n_dofs += psi[0].m_data[g].size ();

    // Create the "Shell" Matrix
    // (the context will be a pointer to this Solver Object)
    ierr = MatCreateShell (PETSC_COMM_WORLD, n_dofs, n_dofs, n_dofs, n_dofs,
        this, &shell_mat);
    ierr = MatShellSetOperation (shell_mat, MATOP_MULT, (void
    (*) ()) shell_vmult<dim>);;

    ierr = EPSSetOperators (eps, shell_mat, NULL);
    ierr = EPSSetProblemType (eps, EPS_NHEP); // The problem is non-symmetric
    ierr = EPSSetTolerances (eps, tol_eps, max_iterations_eps);
    ierr = EPSSetWhichEigenpairs (eps, EPS_LARGEST_MAGNITUDE);
    ierr = EPSSetType (eps, EPSKRYLOVSCHUR);

    // Set Up EPS and Solve
    ierr = EPSSetFromOptions (eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSetUp (eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSolve (eps);
    Assert(ierr == 0, ExcMessage("Error solving SLEPc problem."));

//    PetscViewer viewer;
//    PetscViewerCreate (PETSC_COMM_WORLD, &viewer);
//    PetscViewerSetType (viewer, PETSCVIEWERASCII);
//    PetscViewerFileSetMode (viewer, FILE_MODE_WRITE);
//    PetscViewerFileSetName (viewer, "view.txt");
//    ierr = EPSView (eps, viewer);

    // Resize output vectors
    ierr = VecCreateMPI (PETSC_COMM_WORLD, n_dofs, n_dofs, &tmp_vector);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));

    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      ierr = EPSGetEigenpair (eps, eig, &eigenvalues[eig], PETSC_NULL,
          static_cast<Vec> (tmp_vector), PETSC_NULL);
      Assert(ierr == 0, ExcMessage("Error solving getting the problems."));

      copy_to_sn_vector (tmp_vector, psi[eig]);
    }

    ierr = EPSGetIterationNumber (eps, &n_iterations);

    return ierr;
  }

  /**
   *
   */
  template <int dim>
  void
  EigenSolverSLEPC<dim>::copy_to_petsc_vector (const SnVector & sn_vector,
                                               Vec petsc_vector)
  {
    PetscScalar *avec;
    // PetscErrorCode ierr;
    VecGetArray (petsc_vector, &avec);
    int i = 0;
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int b = 0; b < sn_vector.m_data[g].size (); ++b)
      {
        avec[i] = sn_vector.m_data[g][b];
        i++;
      }
    VecRestoreArray (petsc_vector, &avec);
  }

  /**
   *
   */
  template <int dim>
  void
  EigenSolverSLEPC<dim>::copy_to_sn_vector (Vec petsc_vector,
                                            SnVector & sn_vector)
  {
    PetscScalar *avec;
    // PetscErrorCode ierr;
    VecGetArray (petsc_vector, &avec);
    int i = 0;
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int b = 0; b < sn_vector.m_data[g].size (); ++b)
      {
        sn_vector.m_data[g][b] = avec[i];
        i++;
      }
    VecRestoreArray (petsc_vector, &avec);
  }

  /**
   * @brief Function defined that multiplies the shell matrix by a vector.
   */
  template <int dim>
  PetscErrorCode
  shell_vmult (Mat shell_mat,
               Vec src,
               Vec dst)
  {
    PetscErrorCode ierr;

    // The context of the shell matrix is a pointer to the EigenSolverSLEPC object
    // so we can access the data of the problem.
    void *ctx;
    ierr = MatShellGetContext (shell_mat, &ctx);
    CHKERRQ(ierr);
    EigenSolverSLEPC<dim> *EPSobject = (EigenSolverSLEPC<dim>*) ctx;
    /* Copy to SnVector */
    EPSobject->copy_to_sn_vector (src, EPSobject->sn_vector_tmp_src);
    /* Apply the operator y = A x. */
    EPSobject->m_eigen_prob->vmult (EPSobject->sn_vector_tmp_dst,
        EPSobject->sn_vector_tmp_src);
    /* Copy back to Petsc Vector */
    EPSobject->copy_to_petsc_vector (EPSobject->sn_vector_tmp_dst, dst);

    EPSobject->n_matvec_multiplications++;
    log_print_text (
        "SLEPc Multiplication: " + num_to_string (
            EPSobject->n_matvec_multiplications), 1);

    CHKERRQ(ierr);
    return 0;
  }

  // Explicit Specializations
  template class EigenSolverSLEPC<1> ;
  template class EigenSolverSLEPC<2> ;
  template class EigenSolverSLEPC<3> ;

  template PetscErrorCode
  shell_vmult<1> (Mat shell_mat,
                  Vec src,
                  Vec dst);
  template PetscErrorCode
  shell_vmult<2> (Mat shell_mat,
                  Vec src,
                  Vec dst);
  template PetscErrorCode
  shell_vmult<3> (Mat shell_mat,
                  Vec src,
                  Vec dst);

}
// end of namespace Forest
