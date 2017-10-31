/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op_matrix_built.h
 * @brief  TransOpMatrixBuilt class template declarations
 */

#ifndef FOREST_TRANS_OP_MATRIX_BUILT_H
#define FOREST_TRANS_OP_MATRIX_BUILT_H

#include "neutronics/trans_op.h"

#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&
#include "neutronics/path_to_intergridmap_dealii.h"

#include "deal.II/base/point.h"
#include "deal.II/lac/precondition.h"
#include "deal.II/lac/precondition_block.h"
#include "deal.II/lac/solver_control.h"
#include "deal.II/lac/solver_gmres.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/lac/vector.h"
#include "deal.II/lac/constraint_matrix.h"  // for ConstraintMatrix
#include "deal.II/lac/sparsity_pattern.h"  // for SparsityPattern
#include "deal.II/meshworker/dof_info.h"   // for DoFInfo
#include "deal.II/meshworker/integration_info.h" // for IntegrationInfo
#include "deal.II/meshworker/local_results.h"         // for DoFHandler

#include <vector>
#include <memory>

#ifndef DOXYGEN
namespace dealii { template <class GridClass> class InterGridMapSebas; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   *
   * @class TransOpMatrixBuilt
   *
   * @ingroup ForestNeutronics
   *
   * @brief Transport operator and sweep, matrix-built version
   *
   * @details @copydetails TransOp
   *
   * @todo Add details about the matrix-built implementation
   */

  template <int dim>
  class TransOpMatrixBuilt : public TransOp<dim>
  {
  public:

    /**
     * @brief Constructor
     * @param state
     * @param ord
     */
    TransOpMatrixBuilt (State<dim> & state,
                        std::shared_ptr<QuadratureBase<dim> >& ord);

    /** @brief Destructor */
    virtual
    ~TransOpMatrixBuilt ()
    {
    }

    /* Declaration of inherited pure virtual methods */

    virtual
    void
    apply_L0_inv_g_iord (Vector<double> & phi,
                         const Vector<double> & source,
                         unsigned int g,
                         unsigned int i_ord);

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
      return 0.0;
    };
    virtual
    double
    get_sweep_walltime() const
    {
      return 0.0;
    };

  private:

    /** @brief Build all the matrices and preconditioners */
    void
    build ();

    void
    setup_angular_dof_handler ();

    void
    setup_system ();

    void
    build_trans_matrix ();

    void
    build_prec_matrix ();

    void
    apply_L0_inv_matrix_built (Vector<double> & phi,
                               const Vector<double> & source,
                               unsigned int g,
                               unsigned int i_ord);

    // the following typedef and functions is to use the MeshWorker

    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    void
    trans_integrate_cell_term (DoFInfo &dinfo,
                               CellInfo &info,
                               unsigned int g,
                               const Point<dim> & beta,
                               double weight);

    void
    integrate_influx_boundary_term (DoFInfo &dinfo,
                                    CellInfo &info,
                                    const Point<dim> & beta,
                                    unsigned int boundary_face);

    void
    integrate_outflux_boundary_term (DoFInfo &dinfo,
                                     CellInfo &info,
                                     const Point<dim> & beta);

    void
    integrate_face_term (DoFInfo &dinfo1,
                         DoFInfo &dinfo2,
                         CellInfo &info1,
                         CellInfo &info2,
                         const Point<dim> &beta);

    /**
     * @brief Project scalar to angular mesh
     *
     * Projecting the mesh from
     * dof_handler to dof_handler_angular, and then
     * we solve the system,  where the system is built
     * with the order in dof_handler_angular
     *
     * @note The transformation has been calculated a priory
     * and then we the function "interpolate_to_different_mesh"
     * with it (so we can save time at the cost of some memory
     * allocation)
     *
     */
    void
    project_scalar_to_angular_mesh (unsigned int i,
                                    const Vector<double> & vec_scalar,
                                    Vector<double> & vec_angular);

    /**
     * @brief Project angular to scalar mesh
     *
     * Projecting the mesh from
     * dof_handler_angular to dof_handler
     * after solving the system
     *
     * @note The transformation has been calculated a priory
     * and then we the function "interpolate_to_different_mesh"
     * with it (so we can save time at the cost of some memory
     * allocation)
     *
     */
    void
    project_angular_to_scalar_mesh (unsigned int i,
                                    const Vector<double> & vec_angular,
                                    Vector<double> & vec_scalar);

    const State<dim> & mp_state;
    std::shared_ptr<QuadratureBase<dim> > & mp_ord;
    const DoFHandler<dim, dim> & mp_dof_handler;

    std::vector<DoFHandler<dim, dim> > m_dof_handler_angular;

    std::vector<std::vector<unsigned int> > m_renumbering_dof_indices;
    std::vector<std::vector<unsigned int> > m_reverse_dof_indices;

    std::vector<SparsityPattern> m_trans_pattern;

    /**
     * @brief Transport matrix
     *
     * This matrix have \c n_groups blocks, each
     * block with dimensions \f$ (mn)\times (mn) \f$,
     * being \c m (the number of discrete ordinates
     * directions) the number of blocks in  \c trans_matrix,
     * and \f$ n \f$ the number of spatial unknowns.
     */

    std::vector<std::vector<SparseMatrix<double> > > m_trans_matrix;

    //std::vector<std::vector<PreconditionBlockSOR<SparseMatrix<double>, double> > > m_prec_matrix;
    std::vector<std::vector<PreconditionBlock<SparseMatrix<double>, double> > > m_prec_matrix;
//    std::vector<std::vector<PreconditionBlockSOR<SparseMatrix<double>, double> > > m_prec_matrix;
//
//    std::vector<std::vector<PreconditionBlockSOR<SparseMatrix<double>, double> > > m_prec_matrix;
//    std::vector<std::vector<SparseILU<double> > > m_prec_matrix;
//    std::vector<std::vector<TrilinosWrappers::PreconditionAMG> > m_prec_matrix;

    std::vector<Vector<double> > m_projected_block_dst;
    std::vector<Vector<double> > m_projected_block_src;

    ConstraintMatrix constraints;
    std::vector<InterGridMapSebas<DoFHandler<dim, dim> > > m_grid_m_to_d_map;
    std::vector<InterGridMapSebas<DoFHandler<dim, dim> > > m_grid_d_to_m_map;

  };

} // end of namespace Forest

#endif /* FOREST_TRANS_OP_MATRIX_BUILT_H */
