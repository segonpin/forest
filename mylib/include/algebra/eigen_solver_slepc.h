/**
 * @author  Antoni Vidal. UPV, 2017.
 * @file    algebra/eigen_solver_slepc.h
 * @date    Modified on 31/01/2017
 * @brief   EigenSolver class to solve the function with SLEPCS
 */

#ifndef FOREST_EIGEN_SOLVER_SLEPC_H
#define FOREST_EIGEN_SOLVER_SLEPC_H

#include "algebra/eigen_prob.h"
#include "algebra/eigen_prob_factory.h"

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/slepc_solver.h>

#include <slepceps.h>
#include <petscksp.h>

#include <memory>

#ifndef DOXYGEN
namespace Forest
{
  template <int dim>
  class ProbGeom;
}
namespace Forest
{
  class Input;
}
#endif

#include "algebra/sn_vector.h"

namespace Forest
{

  /**
   *
   */
  template <int dim>
  PetscErrorCode
  shell_vmult (Mat PescMat,
               Vec src,
               Vec dst);

  /**
   * @class EigenSolver
   * @ingroup ForestAlgebra
   * @brief Eigenvalue problem solver
   *
   * @details
   *
   * A detailed description can be found in "ref" lattice_sec
   */
  template <int dim>
  class EigenSolverSLEPC
  {
  public:

    /**
     * @brief Constructor. It needs the eigenvalues problem pointer
     *  (to get the matrix vector product), the number of eigenvalues required
     *  and the number of column vectors involved in the computation
     *   (Krylov subspace dimension).
     */
    EigenSolverSLEPC (std::shared_ptr<EigenProb<dim> > eigen_prob,
                      const unsigned int n_eigenvalues,
                      const unsigned int n_cv);

    /**
     *
     */
    int
    solve (std::vector<double> &eigenvalues,
           std::vector<SnVector> & psi);

    /**
     *
     */
    void
    copy_to_petsc_vector (const SnVector & sn_vector,
                          Vec petsc_vector);

    /**
     *
     */
    void
    copy_to_sn_vector (Vec petsc_vector,
                       SnVector & sn_vector);

    ~EigenSolverSLEPC ();

    /** @brief Container for the eigenproblem. */
    std::shared_ptr<EigenProb<dim> > m_eigen_prob;

    int n_iterations;
    int n_matvec_multiplications;

    double tol_eps;
    unsigned int max_iterations_eps;


    Vec tmp_vector;
    SnVector sn_vector_tmp_src, sn_vector_tmp_dst;

  private:

    Mat shell_mat;
    EPS eps;

    const unsigned int n_eigenvalues;
    const unsigned int n_cv;
    unsigned int n_groups;
  };

} // end of namespace Forest

#endif /* FOREST_EIGEN_SOLVER_SLEPC_H */
