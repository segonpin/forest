/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op_matrix_built.h
 * @brief  DiffOpMatrixBuilt class template declarations
 */

#ifndef FOREST_DIFF_OP_MATRIX_BUILT_H
#define FOREST_DIFF_OP_MATRIX_BUILT_H

#include "neutronics/diff_op.h"
#include "neutronics/state_fwd.h"       // For class State&
#include "neutronics/matrix_integrator.h"

#include "deal.II/lac/block_sparse_matrix.h"
#include "deal.II/fe/mapping_q1.h"      // for MappingQ1
#include "deal.II/lac/precondition_block.h"
#include "deal.II/lac/precondition.h"
#include "deal.II/lac/sparse_ilu.h"
#include "deal.II/lac/sparse_mic.h"

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
   * @class DiffOpMatrixBuilt
   *
   * @ingroup ForestNeutronics
   *
   * @brief Diffusion operator and sweep, matrix-built version
   *
   * @details @copydetails DiffOp
   *
   * @todo Add details about the matrix-built implementation
   */

  template <int dim>
  class DiffOpMatrixBuilt : public DiffOp<dim>
  {
  public:

    /**
     * @brief Constructor
     * @param state
     */
    DiffOpMatrixBuilt (State<dim> & state);

    /** @brief Destructor */
    virtual
    ~DiffOpMatrixBuilt ()
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
    vmult_matrix_built_g_h (BlockVector<double> & phi,
                            const BlockVector<double> & source,
                            const unsigned int g);

    const State<dim> & mp_state;
    const FiniteElement<dim, dim> & mp_fe;
    const MappingQ1<dim> & mp_mapping;
    const DoFHandler<dim, dim> & mp_dof_handler;

    BlockSparsityPattern m_Diff_pattern;
    std::vector<BlockSparseMatrix<double> > m_Diff_matrix;


    //std::vector<PreconditionJacobi<SparseMatrix<double> > > m_prec;
    std::vector<PreconditionBlockJacobi<SparseMatrix<double>, double> > m_prec;
    //std::vector<PreconditionBlockSSOR<SparseMatrix<double>, double> > m_prec;
    //std::vector<SparseILU<double> > m_prec;
    //std::vector<SparseMIC<double> > m_prec;

  };

} // end of namespace Forest

#endif /* FOREST_DIFF_OP_MATRIX_BUILT_H */
