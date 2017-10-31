/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_matrix_built.h
 * @brief  ScatteringMatrixBuilt class template declarations
 */

#ifndef FOREST_SCATTERING_MATRIX_BUILT_H
#define FOREST_SCATTERING_MATRIX_BUILT_H

#include "neutronics/state.h"
#include "neutronics/scattering_source.h"

#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/meshworker/dof_info.h"  // for DoFInfo
#include "deal.II/meshworker/integration_info.h"  // for IntegrationInfo
#include "deal.II/fe/mapping_q1.h"      // for MappingQ1
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/lac/block_sparse_matrix.h" // for BlockSparseMatrix

#include <vector>                       // for vector

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
namespace dealii { template <int dim, int spacedim> class FiniteElement; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @brief Scattering Source, matrix-built version
   * @ingroup ForestNeutronics
   */

  template <int dim>
  class ScatteringMatrixBuilt : public ScatteringSource<dim>
  {

  public:

    /** @brief Constructor */
    ScatteringMatrixBuilt (State<dim> & state);

    /** @brief Destructor */
    virtual
    ~ScatteringMatrixBuilt ()
    {
    }

    /* Declaration of inherited pure virtual methods */

    virtual bool is_matrix_free()
    {
      return false;
    };

    virtual
    double
    memory_consumption () const;

    /**
     @brief This is the building block for the functions in the class.
     @details Applied in a group wise level, switch between matrix-free or
     non-matrix-free to calculate the scattering from g to h of phi and add it
     to the value of source.
     @param source
     @param phi
     @param g
     @param h
     */
    virtual
    void
    vmult_add_g_h (BlockVector<double> & source,
                   const BlockVector<double> & phi,
                   const unsigned int g,
                   const unsigned int h);

    virtual bool no_upscatt(const unsigned int g) const;

  private:

    typedef MeshWorker::DoFInfo<dim> DoFInfo;
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    const State<dim> & mp_state;
    const FiniteElement<dim, dim> & mp_fe;
    const MappingQ1<dim> & mp_mapping;
    const DoFHandler<dim, dim> & mp_dof_handler;

    BlockSparsityPattern m_scatt_pattern;
    std::vector<std::vector<BlockSparseMatrix<double> > > m_scatt_matrix;

    void
    vmult_add_g_h_matrix_built (BlockVector<double> & source,
                   const BlockVector<double> & phi,
                   const unsigned int g,
                   const unsigned int h);

    /**
     @brief Prepare sparsity patterns to setup the systems.
     */
    void
    setup_system ();

    /**
     @brief Build the system for the non-matrix-free option.
     */
    void
    build_system ();

    /**
     @brief Integrates over the cell with a given scattering coefficient.
     @details Used to build the scattering matrix.
     @param dinfo
     @param info
     @param g
     @param h
     */
    void
    scatt_integrate_cell_term (DoFInfo &dinfo,
                               CellInfo &info,
                               const unsigned int g,
                               const unsigned int h);

    void set_non_zero_pattern();
    void print_non_zero_pattern();

    std::vector<std::vector<bool> > m_non_zero;

  };

} // end of namespace Forest

#endif /* FOREST_SCATTERING_MATRIX_BUILT_H */
