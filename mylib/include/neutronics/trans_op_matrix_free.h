/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_op_matrix_free.h
 * @brief  TransOpMatrixFree class template declarations
 */

#ifndef FOREST_TRANS_OP_MATRIX_FREE_H
#define FOREST_TRANS_OP_MATRIX_FREE_H

#include "neutronics/trans_op.h"
#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"       // For class State&

#include "deal.II/lac/full_matrix.h"     // for FullMatrix
#include "deal.II/fe/fe_values.h"        // for FEFaceValues, etc

#include <memory>
#include <vector>

#ifndef DOXYGEN
namespace dealii { template <typename > class Vector; }
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class TransOpMatrixFree
   *
   * @ingroup ForestNeutronics
   *
   * @brief Transport operator and sweep, matrix-free version
   *
   * @details @copydetails TransOp
   *
   * @todo Add details about the matrix-free implementation
   */

  template <int dim>
  class TransOpMatrixFree : public TransOp<dim>
  {
  public:

    /**
     * @brief Constructor
     * @param state
     * @param ord
     */
    TransOpMatrixFree (State<dim> & state,
                       std::shared_ptr<QuadratureBase<dim> > & ord);

    /** @brief Destructor */
    virtual
    ~TransOpMatrixFree ()
    {
    }

    /* Declaration of inherited pure virtual methods */

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
      return 0.0;
    };

    virtual
    double
    get_sweep_walltime() const
    {
      return 0.0;
    };

  private:

    void
    apply_L0_inv_matrix_free (Vector<double> & phi,
                              const Vector<double> & source,
                              const unsigned int g,
                              const unsigned int i_ord);

    void
    apply_L0_inv_mf_pre0 (Vector<double> & phi,
                            const Vector<double> & source,
                            const unsigned int g,
                            const unsigned int i_ord);

    void
    apply_L0_inv_mf_pre (Vector<double> & phi,
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

  };

} // end of namespace Forest

#endif /* FOREST_TRANS_OP_MATRIX_FREE_H */
