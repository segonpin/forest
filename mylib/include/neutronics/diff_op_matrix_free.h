/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op_matrix_free.h
 * @brief  DiffOpMatrixFree class template declarations
 */

#ifndef FOREST_DIFF_OP_MATRIX_FREE_H
#define FOREST_DIFF_OP_MATRIX_FREE_H

#include "neutronics/diff_op.h"
#include "neutronics/state_fwd.h"       // For class State&
#include "neutronics/matrix_integrator.h"
#include "neutronics/matrix_free_integrator.h"

#include "deal.II/lac/block_sparse_matrix.h"
#include "deal.II/fe/mapping_q1.h"      // for MappingQ1

#include <memory>
#include <vector>

#ifndef DOXYGEN
namespace dealii { template <typename> class Vector; }
namespace dealii { template <int dim, int spacedim>  class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class DiffOpMatrixFree
   *
   * @ingroup ForestNeutronics
   *
   * @brief Diffusion operator and sweep, matrix-free version
   *
   * @details @copydetails DiffOp
   *
   * @todo Add details about the matrix-free implementation
   */

  template <int dim>
  class DiffOpMatrixFree : public DiffOp<dim>
  {
  public:

    /**
     * @brief Constructor
     * @param state
     */
    DiffOpMatrixFree (State<dim> & state);

    /** @brief Destructor */
    virtual
    ~DiffOpMatrixFree ()
    {
    }

    /* Declaration of inherited pure virtual methods */

    virtual
    void
    vmult (BlockVector<double> & phi,
           const BlockVector<double> & source,
           const unsigned int g);

    virtual
    void
    prec (BlockVector<double> & phi,
          const BlockVector<double> & source,
          const unsigned int g) const;

    virtual
    unsigned int
    get_n_groups () const;

    virtual
    double
    memory_consumption () const;

  private:

    void
    setup_system ();
    void
    assemble_system ();

    void
    setup_prec ();
    void
    assemble_prec ();

    void
    vmult_matrix_free_g_h (BlockVector<double> & phi,
                           const BlockVector<double> & source,
                           const unsigned int g);

    const State<dim> & mp_state;
    const FiniteElement<dim, dim> & mp_fe;
    const MappingQ1<dim> & mp_mapping;
    const DoFHandler<dim, dim> & mp_dof_handler;

    std::vector<std::shared_ptr<MatrixFreeIntegrator<dim> > > matfree_integrator;
    // BlockSparsityPattern m_Diff_pattern;
    // std::vector<BlockSparseMatrix<double> > m_Diff_matrix;
    // std::vector<BlockVector<double> >  m_prec_matrix;

    /** @brief shortcut for the cell_iterator type */
    //typedef typename DoFHandler<dim, dim>::active_cell_iterator c_iter;
  };

} // end of namespace Forest

#endif /* FOREST_DIFF_OP_MATRIX_FREE_H */
