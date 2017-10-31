/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op_matrix_free.h
 * @brief  TransOpMatrixFree class template declarations
 */

#ifndef TRANS_OP_MATRIX_FREE_PRE_H
#define TRANS_OP_MATRIX_FREE_PRE_H

#include "neutronics/trans_op.h"

#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#include "deal.II/lac/full_matrix.h"     // for FullMatrix
#include "deal.II/fe/fe_values.h"        // for FEFaceValues, etc

#include <memory>
#include <vector>

#ifndef DOXYGEN
namespace dealii {template <typename > class Vector; }
namespace dealii {template <int dim, int spacedim> class DoFHandler;}
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @brief QZ solver for the local matrices.
   * @details It performs the factorization, and then use this factorization
   * to solve different systems
   * */
  class QZsolver
  {
  public:
    /** here we asume matD and matE in vector format, column major ordering */
    QZsolver (const int &n,
              std::vector<double> &matD,
              std::vector<double> &matE);
    /** here we asume matD and matE in FullMatrix format, row major ordering */
    QZsolver (const int &n,
              FullMatrix<double> &matD,
              FullMatrix<double> &matE);
    ~QZsolver();

    void
    solve (const double sigma,
           std::vector<double> &sol,
           std::vector<double> &b);
    void
    solve (const double sigma,
           Vector<double> &sol,
           Vector<double> &b);
  private:
    //--------------------------------------------------------------------------
    // local variables
    int n, info;
    std::vector<double> matD;
    std::vector<double> matE;
    double* matZ;
    double* matQ;
    std::vector<double> newMat;
    std::vector<double> tmp;
    std::vector<double> solve_work;
    // temprary vector
    // We do not use this, but we need it because of the functions signature
    int m_sdim;
    double* m_alphar;
    double* m_alphai;
    double* m_beta;
    int* m_bwork;
    //--------------------------------------------------------------------------
    // functions
    int
    dgges (const int& n,
           double* matG,
           double* matM,
           double* matQ,
           double* matZ);
    int
    dlaqtr (const int& n,
            const double* t,
            double &scale,
            double* x,
            double* work);
    void
    dgemv (const int& n,
           const double* A,
           std::vector<double> &x,
           std::vector<double> &y);
    void
    dgemv (const int& n,
           const double* A,
           std::vector<double> &x,
           Vector<double> &y);
    void
    dgemvT (const int& n,
            const double* A,
            std::vector<double> &x,
            std::vector<double> &y);
    void
    dgemvT (const int& n,
            const double* A,
            Vector<double> &x,
            std::vector<double> &y);
  };

  /**
   * @brief QZ solver for the local matrices.
   * @details It performs the factorization, and then use this factorization
   * to solve different systems
   * */
  class LUsolver
  {
  public:
    /** here we asume matD and matE in vector format, column major ordering */

    LUsolver (const int &n,
              FullMatrix<double> &matD_,
              FullMatrix<double> &matE_);
    void
    solve (const double sigma,
           Vector<double> &sol,
           Vector<double> &b);
  private:
    int n, m_info;
    std::vector<int> ipiv;

    std::vector<double> matD;
    std::vector<double> matE;
    std::vector<double> newMat;
    std::vector<double> tmp;

    int
    dgetrf ();

    int
    dgetrs (Vector<double> &x,
            Vector<double> &y);
  };

  /**
   * @class TransOpMatrixFree
   *
   * @ingroup ForestNeutronics
   *
   * @brief Transport operator and sweep
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
  class TransOpMatrixFreePre : public TransOp<dim>
  {
  public:

    /**
     * @brief Constructor for this class
     * @param state
     * @param ord
     */
    TransOpMatrixFreePre (State<dim> & state,
                          std::shared_ptr<QuadratureBase<dim>> & ord);

    /** @brief Destructor. */
    virtual
    ~TransOpMatrixFreePre ()
    {
    }

    /** @todo document me */
    virtual
    void
    apply_L0_inv_g_iord (Vector<double> & phi,
                         const Vector<double> & source,
                         const unsigned int g,
                         const unsigned int i_ord);

    virtual
    unsigned int
    get_n_groups () const;

    virtual
    unsigned int
    get_n_angles () const;

    virtual
    double
    memory_consumption () const;

    virtual
    double
    get_sweep_cputime() const
    {
      return m_total_sweep_cputime;
    };
    virtual
    double
    get_sweep_walltime() const
    {
      return m_total_sweep_walltime;
    };

  private:

    void
    prepare_local_matrices ();

    void
    apply_L0_inv_mf_pre0 (Vector<double> & phi,
                          const Vector<double> & source,
                          const unsigned int g,
                          const unsigned int i_ord);

    /**
     * @note This function must fixed to work properly. I just have it here to
     * not lose the ideas...
     */
    void
    apply_L0_inv_mf_pre1 (Vector<double> & phi,
                          const Vector<double> & source,
                          const unsigned int g,
                          const unsigned int i_ord);

    /**
     * @note This function must be removed from here
     */
    void
    apply_L0_inv_matrix_free (Vector<double> & phi,
                              const Vector<double> & source,
                              const unsigned int g,
                              const unsigned int i_ord);

    void
    cell_term_mf (const FEValues<dim> &fe_v,
                  const Point<dim> & beta,
                  const double &sigma,
                  FullMatrix<double> &local_matrix) const;

    void
    outgoing_face_term_mf (const FEFaceValuesBase<dim> &fe_v,
                           const Point<dim> &beta,
                           FullMatrix<double> &local_matrix) const;

    void
    incoming_face_term_mf (const FEFaceValuesBase<dim> &fe_v,
                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                           const Point<dim> &beta,
                           Vector<double> &cell_src,
                           const Vector<double> &local_matrix) const;

    const State<dim> & mp_state;
    const std::shared_ptr<QuadratureBase<dim> > & mp_ord;
    const DoFHandler<dim, dim> & mp_dof_handler;
    unsigned int m_n_angles;

    /** @brief shortcut for the cell_iterator type */
    typedef typename DoFHandler<dim, dim>::active_cell_iterator c_iter;
    /** @brief Cell order for the sweep along every direction. */
    std::vector<std::vector<c_iter> > m_sweep_order;
    /** @brief In order to know if a particular face is an incoming face. */
    std::vector<std::vector<std::vector<bool> > > m_incoming_face;

    /** matrices to be precalculated for every direction */
    std::vector<FullMatrix<double>> m_mat_mass;
    std::vector<std::vector<FullMatrix<double> > > m_mat_incoming;
    std::vector<FullMatrix<double>> m_mat_streaming;

    std::vector<std::shared_ptr<QZsolver>> local_qz_system;
    std::vector<std::shared_ptr<LUsolver>> local_lu_system;

    //time
    double m_total_sweep_cputime;
    double m_total_sweep_walltime;

  };

} // end of namespace Forest

#endif
