/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_matrix_free.h
 * @brief  ScatteringMatrixFree class template declarations
 */

#ifndef FOREST_SCATTERING_MATRIX_FREE_H
#define FOREST_SCATTERING_MATRIX_FREE_H

#include "neutronics/scattering_source.h"
#include "algebra/sn_vector.h"             // for SnVector
#include "neutronics/state_fwd.h"

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
   * @brief Scattering Source, matrix-free version
   * @ingroup ForestNeutronics
   */

  template <int dim>
  class ScatteringMatrixFree : public ScatteringSource<dim>
  {

  public:

    /** @brief Constructor */
    ScatteringMatrixFree (State<dim> & state);

    /** @brief Destructor */
    virtual
    ~ScatteringMatrixFree ()
    {
    }

    /* Declaration of inherited pure virtual methods */

    virtual bool
    is_matrix_free ()
    {
      return true;
    }

    virtual
    double
    memory_consumption () const;

    /**
     @brief Scattering from group @p g to @p h for @p phi added to @p source.
     @details This is the main function of the matrix free. It works at group
     level, i.e., it gives
     \f{align*}{
     Q_h = Q_h + S_{g \to h} \Phi_g
     \f}
     @note When you want to get just the scattering coming from group g, then
     you must provide the vector \f$ Q_h \f$ set to zero, (\f$ Q_h = 0 \f$).
     Instead, in general you will provide the vector \f$ Q_h \f$ filled with
     the fission source, or with the outer scattering, \n
     \f$ Q_h = \sum_{g'\neq g} S_{g' \to g} \Phi_{g'} + F_{g}\Phi_{g} \f$.
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

    const State<dim> & mp_state;
    const FiniteElement<dim, dim> & mp_fe;
    const MappingQ1<dim> & mp_mapping;
    const DoFHandler<dim, dim> & mp_dof_handler;

    void
    vmult_add_g_h_matrix_free (BlockVector<double> & source,
                   const BlockVector<double> & phi,
                   const unsigned int g,
                   const unsigned int h);



    void set_non_zero_pattern();
    void print_non_zero_pattern();

    std::vector<std::vector<bool> > m_non_zero;


//    // Not working yet
//    template <int fe_degree>
//    void
//    vmult_add_g_h_k (BlockVector<double> & source,
//                     const BlockVector<double> & phi,
//                     const unsigned int g,
//                     const unsigned int h);

//    void
//    set_coeff (const unsigned int g,
//               const unsigned int h);
//
//    double
//    get_coeff (const unsigned int user_index);
//
//    std::vector<double> coeff;

  };

} // end of namespace Forest

#endif /* FOREST_SCATTERING_MATRIX_FREE_H */
