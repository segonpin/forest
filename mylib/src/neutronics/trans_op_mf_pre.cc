/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op_matrix_free.cc
 * @brief  Implementation of class template TransOpMatrixFree
 */

#include "neutronics/trans_op_mf_pre.h"

#include "utils/forest_utils_dealii.h"
#include "geometry/path_to_gridgenerator_dealii.h"
#include "geometry/prob_geom.h"
#include "input/input.h"
#include "input/input_geom.h"
#include "input/input_mat.h"
#include "angle/quadrature_base.h"
#include "neutronics/state.h"
#include "neutronics/downstream_new.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/base/config.h"          // for DEAL_II_VERSION_GTE
#include "deal.II/base/geometry_info.h"  // for GeometryInfo
#include "deal.II/base/point.h"          // for Point
#include "deal.II/base/tensor.h"         // for Tensor
#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/dofs/dof_tools.h"      // for FiniteElement, Mapping
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc
#include "deal.II/lac/vector.h"          // for Vector
#include "deal.II/lac/full_matrix.h"     // for FullMatrix
#include "deal.II/fe/fe_values.h"        // for FEFaceValues, etc
#include "deal.II/base/timer.h"

#include <utility>                       // for pair
#include <vector>                        // for vector
#include <memory>                        // for shared_ptr

//#include <cblas.h>
#include "cblas.h"
// lapack namespace
namespace lapack
{

  extern "C"
  {
    // if we are running a multi-thread program, we must not use multithread
    // also in openblas, as explained here
    // > https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded
    // we can avoid it with the following call
    // > openblas_set_num_threads(1);
    // void openblas_set_num_threads(int num_threads);

    // DGEMV - compute the matrix vector products
    // y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
    void
    dgemv_ (const char* TRANS,
            const int* M,
            const int* N,
            const double* alpha,
            const double* A,
            const int* LDA,
            double* X,
            const int* INCX,
            const double* beta,
            double* Y,
            const int* INCY);

    // DGGES - compute for a pair of N-by-N real nonsymmetric
    // matrices A, B the generalized eigenvalues, the generalized
    // real Schur form (S,T), optionally, the left and/or right matrices
    // of Schur vectors (VSL and VSR)
    // http://www.netlib.org/lapack/explore-3.1.1-html/dgges.f.html
    void
    dgges_ (const char* jobvsl,
            const char* jobvsr,
            const char* sort,
            int
            (*delztg) (const double*,
                       const double*,
                       const double*),
            const int* n,
            double* a,
            const int* lda,
            double* b,
            const int* ldb,
            int* sdim,
            double* alphar,
            double* alphai,
            double* beta,
            double* vsl,
            const int* ldvsl,
            double* vsr,
            const int* ldvsr,
            double* work,
            const int* lwork,
            int* bwork,
            int* info);

    // DLAQTR - solve the real quasi-triangular system
    // op(T) * p = scale*c
    // http://www.netlib.org/lapack/explore-3.1.1-html/dlaqtr.f.html
    void
    dlaqtr_ (const int* ltran,
             const int* lreal,
             const int* n,
             const double* t,
             const int* ldt,
             const double* b,
             const double* w,
             double* scale,
             double* x,
             double* work,
             int* info);

    // DGETRF perform lu factorization
    void
    dgetrf_ (int* dim1,
             int* dim2,
             double* a,
             int* lda,
             int* ipiv,
             int* info);

    // DGETRS solves a system using the lu factorization calculated with DGETRF
    void
    dgetrs_ (char *TRANS,
             int *N,
             int *NRHS,
             double *A,
             int *LDA,
             int *IPIV,
             double *B,
             int *LDB,
             int *INFO);

  }

  // order_eigs function
  int
  no_order_ (const double* /*ar*/,
             const double* /*ai*/,
             const double* /*b*/)
  {
    return 0;
  }
} // namespace lapack

namespace Forest
{

  using namespace dealii;

  //---------------------------------------------------------------------------

  LUsolver::LUsolver (const int &n,
                      FullMatrix<double> &matD_,
                      FullMatrix<double> &matE_)
      : n (n),
        m_info (0),
        ipiv (n),
        newMat (std::vector<double> (n * n)),
        tmp (std::vector<double> (n))
  {
    // we first move the matrices to the right format (row major ordering)
    matD.resize (n * n);
    matE.resize (n * n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        matD[i + n * j] = matD_ (i, j);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        matE[i + n * j] = matE_ (i, j);
  }

  void
  LUsolver::solve (const double sigma,
                   Vector<double> &sol,
                   Vector<double> &b)
  {
    // first we add the matrices and store it in a new matrix system
    // we can think about using lapack also for this
    for (unsigned int i = 0; i < matD.size (); ++i)
      newMat[i] = matD[i] + sigma * matE[i];
    // print_vector(newMat, "newMat");

    // Use the LAPACK function getrf for calculating the LU factorization.
    m_info = this->dgetrf ();
    if (m_info)
      std::cout << "dgetrf info = " << m_info << std::endl;

    // Use the LAPACK function getrf for solving using the LU factorization.
    m_info = this->dgetrs (sol, b);
    if (m_info)
      std::cout << "dgetrs info = " << m_info << std::endl;
  }

  int
  LUsolver::dgetrf ()
  {
    int info;
    // call dgges for D and E
    // perform the lu decomposition of the matrix in place
    lapack::dgetrf_ (&n, &n, &*newMat.begin (), &n, &*ipiv.begin (), &info);
    return info;
  }

  int
  LUsolver::dgetrs (Vector<double> &sol,
                    Vector<double> &b)
  {
    int info;
    // we define some constants
    char trans = 'N';
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    // move the right hand side to sol
    sol = b;
    // solve the system and return the solution inside sol
    lapack::dgetrs_ (&trans, &n, &nrhs, &*newMat.begin (), &lda,
        &*ipiv.begin (), &*sol.begin (), &ldb, &info);
    return info;
  }

  QZsolver::~QZsolver ()
  {
    delete[] matZ;
    delete[] matQ;
    delete[] m_alphar;
    delete[] m_alphai;
    delete[] m_beta;
    delete[] m_bwork;
  }

  QZsolver::QZsolver (const int &n,
                      std::vector<double> &matD_,
                      std::vector<double> &matE_)
      : n (n),
        info (0),
        matD (matD_),
        matE (matE_),
        matZ (new double[n * n]),
        matQ (new double[n * n]),
        newMat (std::vector<double> (n * n)),
        tmp (std::vector<double> (n)),
        solve_work (std::vector<double> (n)),
        m_sdim (0),
        m_alphar (new double[n]),
        m_alphai (new double[n]),
        m_beta (new double[n]),
        m_bwork (new int[n])
  {
    info = this->dgges (n, &*matD.begin (), &*matE.begin (), matQ, matZ);
    //std::cout << "info = " << info << std::endl;
  }

  QZsolver::QZsolver (const int &n,
                      FullMatrix<double> &matD_,
                      FullMatrix<double> &matE_)
      : n (n),
        info (0),
        matZ (new double[n * n]),
        matQ (new double[n * n]),
        newMat (std::vector<double> (n * n)),
        tmp (std::vector<double> (n)),
        solve_work (std::vector<double> (n)),
        m_sdim (0),
        m_alphar (new double[n]),
        m_alphai (new double[n]),
        m_beta (new double[n]),
        m_bwork (new int[n])
  {
    // we first move the matrices to the right format (row major ordering)
    matD.resize (n * n);
    matE.resize (n * n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        matD[i + n * j] = matD_ (i, j);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        matE[i + n * j] = matE_ (i, j);
    // now we perform the qz decomposition
    info = this->dgges (n, &*matD.begin (), &*matE.begin (), matQ, matZ);
    //std::cout << "info = " << info << std::endl;
  }

  void
  QZsolver::solve (const double sigma,
                   std::vector<double> &sol,
                   std::vector<double> &b)
  {
    //std::cout << "Solving for sigma = " << sigma << std::endl;
    // first we add the matrices and store it in a new matrix system
    // we can think about using lapack also for this
    for (unsigned int i = 0; i < matD.size (); ++i)
      newMat[i] = matD[i] + sigma * matE[i];
    // print_vector(newMat, "newMat");
    // tmp = Q'*b;
    this->dgemvT (n, matQ, b, tmp);
    // print_vector(tmp,"tmp");
    // sol = inv(M)*tmp
    double scale;
    //double* work = new double[n];
    info = this->dlaqtr (n, &*newMat.begin (), scale, &*tmp.begin (), &*solve_work.begin());
    // I am not doing anything with the scale. If I do not use it, then
    // I should hide it from the interface.
    //std::cout << "scale is = " << scale << std::endl;
    // print_vector(tmp,"tmp");
    // tmp = Z*b;
    this->dgemv (n, matZ, tmp, sol);
    // print_vector(sol,"sol"); // print
  }

  void
  QZsolver::solve (const double sigma,
                   Vector<double> &sol,
                   Vector<double> &b)
  {
    //std::cout << "Solving for sigma = " << sigma << std::endl;
    // first we add the matrices and store it in a new matrix system
    // we can think about using lapack also for this
    for (unsigned int i = 0; i < matD.size (); ++i)
      newMat[i] = matD[i] + sigma * matE[i];
    // print_vector(newMat, "newMat");
    // tmp = Q'*b;
    this->dgemvT (n, matQ, b, tmp);
    // print_vector(tmp,"tmp");
    // sol = inv(M)*tmp
    double scale;
    //double* work = new double[n];
    info = this->dlaqtr (n, &*newMat.begin (), scale, &*tmp.begin (), &*solve_work.begin());
    // I am not doing anything with the scale. If I do not use it, then
    // I should hide it from the interface.
    //std::cout << "scale is = " << scale << std::endl;
    // print_vector(tmp,"tmp");
    // tmp = Z*b;
    this->dgemv (n, matZ, tmp, sol);
    // print_vector(sol,"sol"); // print
  }

  int
  QZsolver::dgges (const int& n,
                   double* matG,
                   double* matM,
                   double* matQ,
                   double* matZ)
  {
    // call dgges for D and E
    int info;
    // We call this few times, so it is ok to not reuse work vectors
    int lwork = -1;
    double work_query;
    // we use this options
    const char jobvsl = 'V';
    const char jobvsr = 'V';
    const char sort = 'N';
    // First call to determine memory needed by work
    lapack::dgges_ (&jobvsl, &jobvsr, &sort, lapack::no_order_, &n, matG,
        &n, matM, &n, &m_sdim, m_alphar, m_alphai, m_beta, matQ, &n, matZ, &n,
        &work_query, &lwork, m_bwork, &info);
    lwork = (int) work_query;
    double* work = new double[lwork];
    // Second call to actually perform the qz decomposition
    lapack::dgges_ (&jobvsl, &jobvsr, &sort, lapack::no_order_, &n, matG,
        &n, matM, &n, &m_sdim, m_alphar, m_alphai, m_beta, matQ, &n, matZ, &n,
        work, &lwork, m_bwork, &info);
    delete[] work;
    return info;
  }

  int
  QZsolver::dlaqtr (const int& n,
                    const double* t,
                    double &scale,
                    double* x,
                    double* work)
  {
    int info;
    const int ltran = 0;
    const int lreal = 1;
    // we do not use b and w because it is for the complex solver
    double b;
    double w;
    // We solve the system
    lapack::dlaqtr_ (&ltran, &lreal, &n, t, &n, &b, &w, &scale, x, work, &info);
    return info;
  }

  void
  QZsolver::dgemv (const int& n,
                   const double* A,
                   std::vector<double> &x,
                   std::vector<double> &y)
  {
    const char trans = 'N'; // no transpose
    const int incx = 1; // increment for x
    const int incy = 1; // increment for y
    const double alpha = 1.0;
    const double beta = 0.0;
    lapack::dgemv_ (&trans, &n, &n, &alpha, A, &n, &*x.begin (), &incx, &beta,
        &*y.begin (), &incy);
  }

  void
  QZsolver::dgemv (const int& n,
                   const double* A,
                   std::vector<double> &x,
                   Vector<double> &y)
  {
    const char trans = 'N'; // no transpose
    const int incx = 1; // increment for x
    const int incy = 1; // increment for y
    const double alpha = 1.0;
    const double beta = 0.0;
    lapack::dgemv_ (&trans, &n, &n, &alpha, A, &n, &*x.begin (), &incx, &beta,
        &*y.begin (), &incy);
  }

  void
  QZsolver::dgemvT (const int& n,
                    const double* A,
                    std::vector<double> &x,
                    std::vector<double> &y)
  {
    const char trans = 'T'; // no transpose
    const int incx = 1; // increment for x
    const int incy = 1; // increment for y
    const double alpha = 1.0;
    const double beta = 0.0;
    lapack::dgemv_ (&trans, &n, &n, &alpha, A, &n, &*x.begin (), &incx, &beta,
        &*y.begin (), &incy);
  }

  void
  QZsolver::dgemvT (const int& n,
                    const double* A,
                    Vector<double> &x,
                    std::vector<double> &y)
  {
    const char trans = 'T'; // no transpose
    const int incx = 1; // increment for x
    const int incy = 1; // increment for y
    const double alpha = 1.0;
    const double beta = 0.0;
    lapack::dgemv_ (&trans, &n, &n, &alpha, A, &n, &*x.begin (), &incx, &beta,
        &*y.begin (), &incy);
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  template <int dim>
  TransOpMatrixFreePre<dim>::TransOpMatrixFreePre (State<dim> & state,
                                                   std::shared_ptr<
                                                       QuadratureBase<dim>> & ord)
      : mp_state (state),
        mp_ord (ord),
        mp_dof_handler (mp_state.get_dof_handler ()),
        m_total_sweep_cputime(0.0),
        m_total_sweep_walltime(0.0)
  {
    deallog << "We have enough mesh regularity to use QZ optimization" << std::endl;
    //--------------------------------------------------------------------------
    // Update the value of n_angles to deal with transport.
    m_n_angles = mp_ord->get_n_angles ();
    const unsigned int n_cells = mp_state.mp_geom.m_tria.n_active_cells ();
    //--------------------------------------------------------------------------
    // We set the sweep order using the particular mesh and each direction
    // and flag the incoming boundaries for every cell at every direction.
    m_sweep_order.resize (m_n_angles, std::vector<c_iter> (n_cells));
    m_incoming_face.resize (m_n_angles,
        std::vector<std::vector<bool> > (n_cells,
            std::vector<bool> (GeometryInfo<dim>::faces_per_cell, false)));
    for (unsigned int n = 0; n < m_n_angles; ++n)
    {
      const Point<dim> beta = vector_to_point<dim> (ord->get_q (n));
      DoFRenumbering::compute_downstream_cells (m_sweep_order[n],
          m_incoming_face[n], mp_dof_handler, beta);
    }
    //--------------------------------------------------------------------------
    // We initialize the local matrices for the different angles
    prepare_local_matrices ();
    // constant
    const unsigned int dofs_per_cell = mp_state.get_fe ().dofs_per_cell;
    // qz
    for (unsigned int n = 0; n < m_n_angles; ++n)
      local_qz_system.push_back (
          std::make_shared<QZsolver> (dofs_per_cell, m_mat_streaming[n],
              m_mat_mass[n]));
    // lu
    for (unsigned int n = 0; n < m_n_angles; ++n)
      local_lu_system.push_back (
          std::make_shared<LUsolver> (dofs_per_cell, m_mat_streaming[n],
              m_mat_mass[n]));
    //--------------------------------------------------------------------------
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::prepare_local_matrices ()
  {
    // short names
    typedef typename DoFHandler<dim, dim>::active_cell_iterator ac_iter;
    typedef typename DoFHandler<dim, dim>::cell_iterator c_iter;
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;
    // we get the finite elements and the degree
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const Mapping<dim> & mapping (mp_state.get_mapping ());
    // we get some constants here
    const unsigned int fe_degree = fe.degree;
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    //-------------------------------------------------------------------------
    // Reinit the following matrices
    // mass matrix (should be allways the same)
    m_mat_mass.resize (m_n_angles,
        FullMatrix<double> (dofs_per_cell, dofs_per_cell));
    // matrices for the edges (this are mass matrices)
    m_mat_incoming.resize (m_n_angles,
        std::vector<FullMatrix<double>> (GeometryInfo<dim>::faces_per_cell,
            FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
    // streaming matrix
    m_mat_streaming.resize (m_n_angles,
        FullMatrix<double> (dofs_per_cell, dofs_per_cell));
    //-------------------------------------------------------------------------
    // initialize the Quadratures, FEValues, FEFaceVAlues and FESubFaceValues
    QGauss<dim> quadrature (fe_degree + 1);
    const UpdateFlags update_flags = update_values | update_gradients
                                     | update_JxW_values
                                     | update_inverse_jacobians
                                     | update_jacobians
                                     | update_quadrature_points;
    FEValues<dim> fe_v (mapping, fe, quadrature, update_flags);
    QGauss<dim - 1> face_quadrature (fe_degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
                                          | update_inverse_jacobians
                                          | update_jacobians
                                          | update_quadrature_points
                                          | update_normal_vectors;
    FEFaceValues<dim> fe_v_face (mapping, fe, face_quadrature,
        face_update_flags);
    const UpdateFlags neighbor_update_flags = update_values;
    FEFaceValues<dim> fe_v_face_neighbor (mapping, fe, face_quadrature,
        neighbor_update_flags);
    //-------------------------------------------------------------------------
    // FEValues for reference cells and faces, with their initialization
    FEValues<dim> fe_v_ref (fe, quadrature,
        update_values | update_gradients | update_JxW_values);
    Triangulation<dim> reference_cell;
    GridGenerator::hyper_cube (reference_cell, 0., 1.);
    fe_v_ref.reinit (reference_cell.begin ());
    // I use shared_ptr because I cannot use FEFaceValues with empty constructor in a vector
    std::vector<std::shared_ptr<FEFaceValues<dim>>>fe_v_f_ref;
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      fe_v_f_ref.push_back (
          std::make_shared<FEFaceValues<dim>> (fe, face_quadrature,
              update_values | update_JxW_values | update_quadrature_points
              | update_normal_vectors));
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      fe_v_f_ref[f]->reinit (reference_cell.begin (), f);
    //-------------------------------------------------------------------------
    // initialize the reference FEValues for a cell surrounded by neighbours
    Triangulation<dim> ref_cells;
    GridGenerator::subdivided_hyper_cube (ref_cells, 3, -1., 2.);
    //const Point<dim> center = vector_to_point<dim>(std::vector<double>(dim,0.5));
    const Point<dim> center =
        (dim == 1) ?
            Point<dim> (0.5) :
            ((dim == 2) ? Point<dim> (0.5, 0.5) : Point<dim> (0.5, 0.5, 0.5));
    DoFHandler<dim> dofs_ref (ref_cells);
    dofs_ref.distribute_dofs (fe);
    ac_iter central_cell;
    for (ac_iter c = dofs_ref.begin_active (); c != dofs_ref.end (); ++c)
      if (c->center ().distance (center) < 1.e-6)
        central_cell = c;
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      fe_v_face.reinit (central_cell, f);
    //-------------------------------------------------------------------------
    // loop over the directions
    for (unsigned int n = 0; n < m_n_angles; ++n)
    {
      // updating the streaming direction
      const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (n));
      m_mat_streaming[n] = 0.0;
      {
        ac_iter cell = m_sweep_order[n][0];
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
          // reset the matrices for the incoming faces before assign new values
          m_mat_incoming[n][f] = 0;
          if (m_incoming_face[n][cell->user_index ()][f])
          {
            fe_v_face.reinit (cell, f);
            // case d), cells are at the same level.
            ac_iter neighbor = central_cell->neighbor (f);
            const unsigned int neighbor2 = central_cell->neighbor_of_neighbor (
                f);
            fe_v_face_neighbor.reinit (neighbor, neighbor2);
            const double beta_n = beta * fe_v_face.boundary_form (0);
            // We use this loop
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  m_mat_incoming[n][f] (i, j) -= beta_n
                      * fe_v_face_neighbor.shape_value (j, q)
                      //* fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);
          }
          else if (not m_incoming_face[n][cell->user_index ()][f])
          {
            // case e), this is an outgoing face.
            fe_v_face.reinit (cell, f);
            const double beta_n = beta * fe_v_face.boundary_form (0);
            // We use this loop
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  m_mat_streaming[n] (i, j) += beta_n
                      * fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);
          }
        }
        fe_v.reinit (cell);
        double detJ = fe_v.jacobian (0).determinant ();
        Tensor<1, dim> betaJ = apply_transformation (
            fe_v.inverse_jacobian (0).transpose (), beta);
        // precalculated streaming matrix in the cell element
        for (unsigned int q = 0; q < fe_v_ref.n_quadrature_points; ++q)
          for (unsigned int i = 0; i < fe_v_ref.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe_v_ref.dofs_per_cell; ++j)
              for (unsigned int d = 0; d < dim; ++d)
                m_mat_streaming[n] (i, j) += -detJ * betaJ[d] * fe_v_ref.JxW (q)
                                             * fe_v_ref.shape_grad (i, q)[d]
                                             * fe_v_ref.shape_value (j, q);
        // precalculated mass matrix in the cell element
        m_mat_mass[n] = 0.0;
        for (unsigned int q = 0; q < fe_v_ref.n_quadrature_points; ++q)
          for (unsigned int i = 0; i < fe_v_ref.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe_v_ref.dofs_per_cell; ++j)
              m_mat_mass[n] (i, j) += detJ * fe_v_ref.JxW (q)
                                      * fe_v_ref.shape_value (j, q)
                                      * fe_v_ref.shape_value (i, q);
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::apply_L0_inv_g_iord (Vector<double> & phi,
                                                  const Vector<double> & source,
                                                  const unsigned int g,
                                                  const unsigned int i_ord)
  {
    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart (); // We restart the timer

    //apply_L0_inv_matrix_free (phi, source, g, i_ord);
    apply_L0_inv_mf_pre0 (phi, source, g, i_ord);
    //apply_L0_inv_mf_pre (phi, source, g, i_ord);

    // Stop timers
    timer.stop ();
    m_total_sweep_cputime+= timer();
    m_total_sweep_walltime+= timer.wall_time ();
  }

  /**
   * @brief Here 0 refers to the number of hanging nodes allowed.
   * It means, no hanging nodes are allowed.
   */
  template <int dim>
  void
  TransOpMatrixFreePre<dim>::apply_L0_inv_mf_pre0 (Vector<double> & dst,
                                                   const Vector<double> & src,
                                                   const unsigned int g,
                                                   const unsigned int i_ord)
  {
    //-------------------------------------------------------------------------
    // short names
    typedef typename DoFHandler<dim, dim>::active_cell_iterator ac_iter;
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;
    // Initialize the destiny vector to zero
    dst = 0;
    // updating the streaming direction is not necessary any more here
    //const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));
    // we get some constants here
    const unsigned int dofs_per_cell = mp_state.get_fe ().dofs_per_cell;
    //-------------------------------------------------------------------------
    // vectors to perform the local operations
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell),
        neighbor_src (dofs_per_cell);
    // to map the dofs from local to global and viceversa
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<unsigned int> neighbor_dof_indices (dofs_per_cell);
    // the next matrices are going to depend on the previous ones
    // and will change depending on the cross section (g,e) and direction (w)
    FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> aux_matrix (dofs_per_cell, dofs_per_cell);
    //-------------------------------------------------------------------------
    for (unsigned int c = 0; c < m_sweep_order[i_ord].size (); ++c)
    {
      // set the cell we are working with.
      ac_iter cell = m_sweep_order[i_ord][c];
      // now we get the source for this cell
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        cell_src (i) = src[local_dof_indices[i]];
      // adding the contributions from the previous cell
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (m_incoming_face[i_ord][cell->user_index ()][f])
        {
          // Case a) If incoming boundary conditions. They have already
          // been applied, otherwise, something as follows must be used
          // fe_v_face.reinit (cell, face_no);
          // boundary_term_matrix_free(fe_v_face, cell_src);
          if (cell->face (f)->at_boundary ())
            continue;
          // case d), cells are at the same level.
          cell->neighbor (f)->get_dof_indices (neighbor_dof_indices);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            neighbor_src (i) = dst[neighbor_dof_indices[i]];
          m_mat_incoming[i_ord][f].vmult_add (cell_src, neighbor_src);
        }
      }
      const bool use_qz = true;
      const bool use_lu = true;
      // solver
      if (use_qz)
      {
        const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
            cell->material_id ()).m_sigmat[g];
        local_qz_system[i_ord]->solve (sigmat, cell_dst, cell_src);
      }
      else if (use_lu)
      {
        const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
            cell->material_id ()).m_sigmat[g];
        local_lu_system[i_ord]->solve (sigmat, cell_dst, cell_src);
      }
      else
      {
        // Construction of the local matrix from precomputed blocks
        local_matrix = 0;
        const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
            cell->material_id ()).m_sigmat[g];
        double sigma = sigmat;
        local_matrix.add (1.0, m_mat_streaming[i_ord], sigma,
            m_mat_mass[i_ord]);
        // we solve the local system
        aux_matrix.invert (local_matrix);
        aux_matrix.vmult (cell_dst, cell_src);
      }

      // we map local-global to put the local solution in the global vector
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        dst[local_dof_indices[i]] += cell_dst (i);
    }
    //-------------------------------------------------------------------------
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::apply_L0_inv_mf_pre1 (Vector<double> & dst,
                                                   const Vector<double> & src,
                                                   const unsigned int g,
                                                   const unsigned int i_ord)
  {
    // short names
    typedef typename DoFHandler<dim, dim>::active_cell_iterator ac_iter;
    typedef typename DoFHandler<dim, dim>::cell_iterator c_iter;
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;

    // Initialize the destiny vector to zero
    dst = 0;

    // updating the streaming direction
    // updating the weight
    const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));

    // we get the finite elements and the degree
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const Mapping<dim> & mapping (mp_state.get_mapping ());
    const unsigned int fe_degree = fe.degree;

    // initialize the Quadratures, FEValues, FEFaceVAlues and FESubFaceValues
    QGauss<dim> quadrature (fe_degree + 1);
    const UpdateFlags update_flags = update_values | update_gradients
                                     | update_JxW_values
                                     | update_inverse_jacobians
                                     | update_jacobians
                                     | update_quadrature_points;
    FEValues<dim> fe_v (mapping, fe, quadrature, update_flags);
    QGauss<dim - 1> face_quadrature (fe_degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
                                          | update_inverse_jacobians
                                          | update_jacobians
                                          | update_quadrature_points
                                          | update_normal_vectors;
    FEFaceValues<dim> fe_v_face (mapping, fe, face_quadrature,
        face_update_flags);
    FESubfaceValues<dim> fe_v_subface (mapping, fe, face_quadrature,
        face_update_flags);
    const UpdateFlags neighbor_update_flags = update_values;
    FEFaceValues<dim> fe_v_face_neighbor (mapping, fe, face_quadrature,
        neighbor_update_flags);

    //------------------------------------------
    FEValues<dim> fe_v_ref (fe, quadrature,
        update_values | update_gradients | update_JxW_values);
    Triangulation<dim> reference_cell;
    GridGenerator::hyper_cube (reference_cell, 0., 1.);
    fe_v_ref.reinit (reference_cell.begin ());

    /*std::vector<FEFaceValues<dim>> fe_v_f_ref (
     GeometryInfo<dim>::faces_per_cell,
     FEFaceValues<dim> (fe, face_quadrature,
     update_values | update_JxW_values | update_quadrature_points
     | update_normal_vectors));
     for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
     fe_v_f_ref[f].reinit (reference_cell.begin (), f);*/

    // I use shared_ptr because I cannot use FEFaceValues with empty constructor in a vector
    std::vector<std::shared_ptr<FEFaceValues<dim>>>fe_v_f_ref;
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      fe_v_f_ref.push_back (
          std::make_shared<FEFaceValues<dim>> (fe, face_quadrature,
              update_values | update_JxW_values | update_quadrature_points
              | update_normal_vectors));
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      fe_v_f_ref[f]->reinit (reference_cell.begin (), f);

    // I use shared_ptr because I cannot use FEFaceValues with empty constructor in a vector
    std::vector<std::shared_ptr<FESubfaceValues<dim>>>fe_v_sub_f_ref;
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      for (unsigned int fsub = 0;
          fsub < GeometryInfo<dim>::max_children_per_face; ++fsub)
        fe_v_sub_f_ref.push_back (
            std::make_shared<FESubfaceValues<dim>> (fe, face_quadrature,
                update_values | update_JxW_values | update_quadrature_points
                | update_normal_vectors));
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      for (unsigned int fsub = 0;
          fsub < GeometryInfo<dim>::max_children_per_face; ++fsub)
        fe_v_sub_f_ref[f * GeometryInfo<dim>::max_children_per_face + fsub]->reinit (
            reference_cell.begin (), f, fsub);

    //------------------------------------------
    // Create local matrices
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> mat_mass (dofs_per_cell, dofs_per_cell);
    std::vector<FullMatrix<double>> mat_incoming (
        GeometryInfo<dim>::faces_per_cell,
        FullMatrix<double> (dofs_per_cell, dofs_per_cell));
    std::vector<FullMatrix<double>> mat_incoming_sub (
        GeometryInfo<dim>::faces_per_cell * GeometryInfo<dim>::max_children_per_face,
        FullMatrix<double> (dofs_per_cell, dofs_per_cell));
    FullMatrix<double> mat_streaming (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> aux_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell),
        neighbor_src (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<unsigned int> neighbor_dof_indices (dofs_per_cell);

    //-------------------------------------------------------------------------
    {
      Triangulation<dim> ref_cell_childs;
      GridGenerator::subdivided_hyper_cube (ref_cell_childs, 3, -1., 2.);
      const Point<dim> center =
          (dim == 1) ?
              Point<dim> (0.5) :
              ((dim == 2) ? Point<dim> (0.5, 0.5) : Point<dim> (0.5, 0.5, 0.5));
      {
        typename Triangulation<dim>::active_cell_iterator cell;
        for (cell = ref_cell_childs.begin_active ();
            cell != ref_cell_childs.end (); ++cell)
          if (cell->center ().distance (center) > 1.e-6)
            cell->set_refine_flag ();
        ref_cell_childs.execute_coarsening_and_refinement ();
      }
      DoFHandler<dim> dofs_ref_child (ref_cell_childs);
      dofs_ref_child.distribute_dofs (fe);

      ac_iter central_cell;
      for (ac_iter cell = dofs_ref_child.begin_active ();
          cell != dofs_ref_child.end (); ++cell)
        if (cell->center ().distance (center) < 1.e-6)
          central_cell = cell;

      const int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

      for (unsigned int f_no = 0; f_no < faces_per_cell; ++f_no)
        if (central_cell->face (f_no)->has_children ())
          for (unsigned int subf_no = 0;
              subf_no < central_cell->face (f_no)->number_of_children ();
              ++subf_no)
            fe_v_subface.reinit (central_cell, f_no, subf_no);
      //fe_v_subface[f_no][subf_no].reinit (central_cell,f_no,subf_no);

      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        ac_iter cell = central_cell;
        f_iter face = cell->face (f);

        if (face->has_children ())
        {
          // case b), face has children.
          const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
          for (unsigned int subf = 0; subf < face->number_of_children ();
              ++subf)
          {
            c_iter neigh_child = cell->neighbor_child_on_subface (f, subf);
            neigh_child->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-----------------------------
            fe_v_subface.reinit (cell, f, subf);
            fe_v_face_neighbor.reinit (neigh_child, neighbor2);
            /*incoming_face_term_mf (fe_v_subface, fe_v_face_neighbor, beta,
             cell_src, neighbor_src);*/

            const double beta_n = beta * fe_v_subface.boundary_form (0);
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  mat_incoming_sub[f] (i, j) -= beta_n
                      * fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);
          }
        }
        else
        {
          std::cout << "warning: face should always have a children"
                    << std::endl;
        }

      }
    }

    mat_streaming = 0.0;
    {
      ac_iter cell = m_sweep_order[i_ord][0];
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        mat_incoming[f] = 0;
        unsigned int fsub = 0;
        mat_incoming_sub[f * GeometryInfo<dim>::max_children_per_face + fsub] =
            0;

        if (m_incoming_face[i_ord][cell->user_index ()][f])
        {
          f_iter face = cell->face (f);

          if (face->has_children ())
          {
            // case b), face has children.
            const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
            for (unsigned int subf = 0; subf < face->number_of_children ();
                ++subf)
            {
              c_iter neigh_child = cell->neighbor_child_on_subface (f, subf);
              neigh_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                neighbor_src (i) = dst[neighbor_dof_indices[i]];
              // Change for a precalculated matrix?-----------------------------
              fe_v_subface.reinit (cell, f, subf);
              fe_v_face_neighbor.reinit (neigh_child, neighbor2);
              /*incoming_face_term_mf (fe_v_subface, fe_v_face_neighbor, beta,
               cell_src, neighbor_src);*/

              const double beta_n = beta * fe_v_subface.boundary_form (0);
              for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                  ++q)
                for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell;
                      ++j)
                    mat_incoming_sub[f] (i, j) -= beta_n
                        * fe_v_f_ref[f]->shape_value (j, q)
                        * fe_v_f_ref[f]->shape_value (i, q)
                        * fe_v_f_ref[f]->JxW (q);
            }
          }
          else if (cell->neighbor_is_coarser (f))
          {
            // case c), neighbor is coarser.
            ac_iter neighbor = cell->neighbor (f);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-------------------------------
            std::pair<unsigned int, unsigned int> neighbor_face_subface =
                cell->neighbor_of_coarser_neighbor (f);
            fe_v_face.reinit (cell, f);
            fe_v_subface.reinit (neighbor, neighbor_face_subface.first,
                neighbor_face_subface.second);
            //incoming_face_term_mf (fe_v_face, fe_v_subface, beta, cell_src,
            //    neighbor_src);

            const double beta_n = beta * fe_v_face.boundary_form (0);
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  mat_incoming_sub[f] (i, j) -= beta_n
                      * fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);

          }
          else
          {
            // case d), cells are at the same level.
            // Precalculated matrix?-----------------------------
            fe_v_face.reinit (cell, f);
            const double beta_n = beta * fe_v_face.boundary_form (0);
            // We use this loop
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  mat_incoming[f] (i, j) -= beta_n
                      * fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);
          }
        }

        if (not m_incoming_face[i_ord][cell->user_index ()][f])
        {
          {
            // case e), this is an outgoing face.
            fe_v_face.reinit (cell, f);

            // @todo Fe_v_face is returning the jacobian of the full cell and not
            // just for the edge. So for the moment, assuming that the cell is the
            // same length along all the axis (then J(d) = J^d) I correct the
            // exponent of the jacobian to be one dimension less
            const double beta_n = beta * fe_v_face.boundary_form (0);
            // We use this loop
            for (unsigned int q = 0; q < fe_v_f_ref[f]->n_quadrature_points;
                ++q)
              for (unsigned int i = 0; i < fe_v_f_ref[f]->dofs_per_cell; ++i)
                for (unsigned int j = 0; j < fe_v_f_ref[f]->dofs_per_cell; ++j)
                  mat_streaming (i, j) += beta_n
                      * fe_v_f_ref[f]->shape_value (j, q)
                      * fe_v_f_ref[f]->shape_value (i, q)
                      * fe_v_f_ref[f]->JxW (q);
          }
        }
      }
      fe_v.reinit (cell);
      double detJ = fe_v.jacobian (0).determinant ();
      Tensor<1, dim> betaJ = apply_transformation (
          fe_v.inverse_jacobian (0).transpose (), beta);
      // precalculated streaming matrix in the cell element
      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < fe_v_ref.dofs_per_cell; ++i)
          for (unsigned int j = 0; j < fe_v_ref.dofs_per_cell; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              mat_streaming (i, j) += -detJ * betaJ[d] * fe_v_ref.JxW (q)
                                      * fe_v_ref.shape_grad (i, q)[d]
                                      * fe_v_ref.shape_value (j, q);
      // precalculated mass matrix in the cell element
      mat_mass = 0.0;
      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < fe_v_ref.dofs_per_cell; ++i)
          for (unsigned int j = 0; j < fe_v_ref.dofs_per_cell; ++j)
            mat_mass (i, j) += detJ * fe_v_ref.JxW (q)
                               * fe_v_ref.shape_value (j, q)
                               * fe_v_ref.shape_value (i, q);
    }

    for (unsigned int cell_i = 0; cell_i < m_sweep_order[i_ord].size ();
        ++cell_i)
    {
      // set the cell we are working with.
      ac_iter cell = m_sweep_order[i_ord][cell_i];

      // Initialize the local matrix to be zero
      local_matrix = 0;
      // Copy the src for the cell
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        cell_src (i) = src[local_dof_indices[i]];
      }

      // The first is to run over all the neighbor cells and apply them
      // when the flux is incoming to this cell
      // the streaming ordering takes care that the fictitious matrix is
      // lower triangular
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {

        if (m_incoming_face[i_ord][cell->user_index ()][f])
        {
          f_iter face = cell->face (f);

          // Case a) If we are at the zero incoming boundary conditions have already
          // been applied, otherwise, something as follows should be used
          // fe_v_face.reinit (cell, face_no);
          // boundary_term_matrix_free(fe_v_face, cell_src);
          if (face->at_boundary ())
          {
            continue;
          }
          else if (face->has_children ())
          {
            // case b), face has children.
            const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
            for (unsigned int subf = 0; subf < face->number_of_children ();
                ++subf)
            {
              c_iter neigh_child = cell->neighbor_child_on_subface (f, subf);
              neigh_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                neighbor_src (i) = dst[neighbor_dof_indices[i]];
              // Change for a precalculated matrix?-----------------------------
              fe_v_subface.reinit (cell, f, subf);
              fe_v_face_neighbor.reinit (neigh_child, neighbor2);
              incoming_face_term_mf (fe_v_subface, fe_v_face_neighbor, beta,
                  cell_src, neighbor_src);
            }
          }
          else if (cell->neighbor_is_coarser (f))
          {
            // case c), neighbor is coarser.
            ac_iter neighbor = cell->neighbor (f);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-------------------------------
            std::pair<unsigned int, unsigned int> neighbor_face_subface =
                cell->neighbor_of_coarser_neighbor (f);
            fe_v_face.reinit (cell, f);
            fe_v_subface.reinit (neighbor, neighbor_face_subface.first,
                neighbor_face_subface.second);
            incoming_face_term_mf (fe_v_face, fe_v_subface, beta, cell_src,
                neighbor_src);
          }
          else
          {
            // case d), cells are at the same level.
            ac_iter neighbor = cell->neighbor (f);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-----------------------------
            const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
            fe_v_face.reinit (cell, f);
            fe_v_face_neighbor.reinit (neighbor, neighbor2);
            incoming_face_term_mf (fe_v_face, fe_v_face_neighbor, beta,
                cell_src, neighbor_src);
          }
        }
        else
        {
          // case e), this is an outgoing face.
          fe_v_face.reinit (cell, f);
          outgoing_face_term_mf (fe_v_face, beta, local_matrix);
        }
      }
      // Now we should build and invert (or solve a system) with the
      // streaming local matrix (with our without the scattering?)
      {
        fe_v.reinit (cell);
        const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
            cell->material_id ()).m_sigmat[g];
        cell_term_mf (fe_v, beta, sigmat, local_matrix);
      }

      /*std::cout << "on  matrix " << std::endl;
       local_matrix.print_formatted (std::cout, 3, true, 9, "0.000e+00");
       std::cout << std::endl << std::endl;*/

      // Construction of the local matrix from precomputed blocks
      local_matrix = 0;
      const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
          cell->material_id ()).m_sigmat[g];
      double sigma = sigmat;
      local_matrix.add (1.0, mat_streaming, sigma, mat_mass);

      /*std::cout << "full matrix " << std::endl;
       local_matrix.print_formatted (std::cout, 3, true, 9, "0.000e+00");
       std::cout << std::endl;
       std::cout << "----------------------------------------" << std::endl;*/

      // we solve/multiply by the block in the diagonal
      // local_matrix.vmult (cell_dst, cell_src);
      aux_matrix.invert (local_matrix);
      aux_matrix.vmult (cell_dst, cell_src);

      // we map from local to global
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        dst[local_dof_indices[i]] += cell_dst (i);
      }
    }
    //std::exit (1);
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::apply_L0_inv_matrix_free (Vector<double> & dst,
                                                       const Vector<double> & src,
                                                       const unsigned int g,
                                                       const unsigned int i_ord)
  {
    // short names
    typedef typename DoFHandler<dim, dim>::active_cell_iterator ac_iter;
    typedef typename DoFHandler<dim, dim>::cell_iterator c_iter;
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;

    // Initialize the destiny vector to zero
    dst = 0;

    // updating the streaming direction
    // updating the weight
    const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));

    // we get the finite elements and the degree
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const Mapping<dim> & mapping (mp_state.get_mapping ());
    const unsigned int fe_degree = fe.degree;

    // initialize the Quadratures, FEValues, FEFaceVAlues and FESubFaceValues
    QGauss<dim> quadrature (fe_degree + 1);
    const UpdateFlags update_flags = update_values | update_gradients
                                     | update_JxW_values
                                     | update_quadrature_points;
    FEValues<dim> fe_v (mapping, fe, quadrature, update_flags);
    QGauss<dim - 1> face_quadrature (fe_degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
                                          | update_quadrature_points
                                          | update_normal_vectors;
    FEFaceValues<dim> fe_v_face (mapping, fe, face_quadrature,
        face_update_flags);
    FESubfaceValues<dim> fe_v_subface (mapping, fe, face_quadrature,
        face_update_flags);
    const UpdateFlags neighbor_update_flags = update_values;
    FEFaceValues<dim> fe_v_face_neighbor (mapping, fe, face_quadrature,
        neighbor_update_flags);

    // Create local matrices
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> aux_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell),
        neighbor_src (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<unsigned int> neighbor_dof_indices (dofs_per_cell);

    for (unsigned int cell_i = 0; cell_i < m_sweep_order[i_ord].size ();
        ++cell_i)
    {
      // set the cell we are working with.
      ac_iter cell = m_sweep_order[i_ord][cell_i];

      // Initialize the local matrix to be zero
      local_matrix = 0;
      // Copy the src for the cell
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        cell_src (i) = src[local_dof_indices[i]];
      }

      // The first is to run over all the neighbor cells and apply them
      // when the flux is incoming to this cell
      // the streaming ordering takes care that the fictitious matrix is
      // lower triangular
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {

        if (m_incoming_face[i_ord][cell->user_index ()][f])
        {
          f_iter face = cell->face (f);

          // Case a) If we are at the zero incoming boundary conditions have already
          // been applied, otherwise, something as follows should be used
          // fe_v_face.reinit (cell, face_no);
          // boundary_term_matrix_free(fe_v_face, cell_src);
          if (face->at_boundary ())
          {
            continue;
          }
          else if (face->has_children ())
          {
            // case b), face has children.
            const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
            for (unsigned int subf = 0; subf < face->number_of_children ();
                ++subf)
            {
              c_iter neigh_child = cell->neighbor_child_on_subface (f, subf);
              neigh_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                neighbor_src (i) = dst[neighbor_dof_indices[i]];
              // Change for a precalculated matrix?-----------------------------
              fe_v_subface.reinit (cell, f, subf);
              fe_v_face_neighbor.reinit (neigh_child, neighbor2);
              incoming_face_term_mf (fe_v_subface, fe_v_face_neighbor, beta,
                  cell_src, neighbor_src);
            }
          }
          else if (cell->neighbor_is_coarser (f))
          {
            // case c), neighbor is coarser.
            ac_iter neighbor = cell->neighbor (f);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-------------------------------
            std::pair<unsigned int, unsigned int> neighbor_face_subface =
                cell->neighbor_of_coarser_neighbor (f);
            fe_v_face.reinit (cell, f);
            fe_v_subface.reinit (neighbor, neighbor_face_subface.first,
                neighbor_face_subface.second);
            incoming_face_term_mf (fe_v_face, fe_v_subface, beta, cell_src,
                neighbor_src);
          }
          else
          {
            // case d), cells are at the same level.
            ac_iter neighbor = cell->neighbor (f);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = dst[neighbor_dof_indices[i]];
            // Change for a precalculated matrix?-----------------------------
            const unsigned int neighbor2 = cell->neighbor_of_neighbor (f);
            fe_v_face.reinit (cell, f);
            fe_v_face_neighbor.reinit (neighbor, neighbor2);
            incoming_face_term_mf (fe_v_face, fe_v_face_neighbor, beta,
                cell_src, neighbor_src);
          }
        }
        else
        {
          // case e), this is an outgoing face.
          fe_v_face.reinit (cell, f);
          outgoing_face_term_mf (fe_v_face, beta, local_matrix);
        }
      }

      // Now we should build and invert (or solve a system) with the
      // streaming local matrix (with our without the scattering?)
      {
        fe_v.reinit (cell);
        const double sigmat = mp_state.mp_data.mp_mat.m_xs.at (
            cell->material_id ()).m_sigmat[g];
        cell_term_mf (fe_v, beta, sigmat, local_matrix);
      }

      // we solve/multiply by the block in the diagonal
      // local_matrix.vmult (cell_dst, cell_src);
      aux_matrix.invert (local_matrix);
      aux_matrix.vmult (cell_dst, cell_src);

      // we map from local to global
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        dst[local_dof_indices[i]] += cell_dst (i);
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::cell_term_mf (const FEValues<dim> &fe_v,
                                           const Point<dim> &beta,
                                           const double &sigma,
                                           FullMatrix<double> &local_matrix) const
  {
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          local_matrix (i, j) += fe_v.shape_value (j, q)
              * (-beta * fe_v.shape_grad (i, q) + sigma
                  * fe_v.shape_value (i, q))
              * JxW[q];
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::outgoing_face_term_mf (const FEFaceValuesBase<dim> &fe_v_face,
                                                    const Point<dim> &beta,
                                                    FullMatrix<double> &local_matrix) const
  {
    const std::vector<Tensor<1, dim> > &normals =
        fe_v_face.get_all_normal_vectors ();
    const std::vector<double> &JxW = fe_v_face.get_JxW_values ();

    // outgoing boundary conditions for the present cell will be
    // added to the local matrix
    for (unsigned int q = 0; q < fe_v_face.n_quadrature_points; ++q)
    {
      const double beta_n = beta * normals[q];
      for (unsigned int i = 0; i < fe_v_face.dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < fe_v_face.dofs_per_cell; ++j)
        {
          local_matrix (i, j) += beta_n * fe_v_face.shape_value (j, q)
                                 * fe_v_face.shape_value (i, q) * JxW[q];
        }
      }
    }
  }

  template <int dim>
  void
  TransOpMatrixFreePre<dim>::incoming_face_term_mf (const FEFaceValuesBase<dim> &fe_v,
                                                    const FEFaceValuesBase<dim> &fe_v_neighbor,
                                                    const Point<dim> &beta,
                                                    Vector<double> &cell_src,
                                                    const Vector<double> &neighbor_src) const
  {
    const std::vector<Tensor<1, dim> > &normals =
        fe_v.get_all_normal_vectors ();
    const std::vector<double> &JxW = fe_v.get_JxW_values ();
    for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
    {
      double neighbor_value = 0;
      for (unsigned int j = 0; j < fe_v_neighbor.dofs_per_cell; ++j)
        neighbor_value += neighbor_src (j) * fe_v_neighbor.shape_value (j, q);
      const double beta_n = beta * normals[q];
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        cell_src (i) -= beta_n * fe_v.shape_value (i, q) * neighbor_value
                        * JxW[q];
    }
  }

  template <int dim>
  unsigned int
  TransOpMatrixFreePre<dim>::get_n_groups () const
  {
    return mp_state.get_n_groups ();
  }

  template <int dim>
  unsigned int
  TransOpMatrixFreePre<dim>::get_n_angles () const
  {
    return m_n_angles;
  }

  template <int dim>
  double
  TransOpMatrixFreePre<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption += 0;
    return memory_consumption;
  }

  template class TransOpMatrixFreePre<1> ;
  template class TransOpMatrixFreePre<2> ;
  template class TransOpMatrixFreePre<3> ;

} // end of namespace Forest
