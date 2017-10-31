/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/moments_to_directions.h
 * @brief   MomentsToDirections class template declarations
 */

#ifndef FOREST_MOMENTS_TO_DIRECTIONS_H
#define FOREST_MOMENTS_TO_DIRECTIONS_H

#include "angle/moment_indexer.h"
#include "angle/quadrature_base_fwd.h"

#include <memory>

#ifndef DOXYGEN
namespace dealii { template <typename> class Vector; }
namespace dealii { template <typename> class BlockVector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @brief Moments to directions class
   * @ingroup ForestAngle
   * @details Given a moment expansion of the angular flux in spherical
   * harmonics, we project it into a discrete ordinates representation
   * of the angular flux,
   *
   * \f{align*}{
   * \Psi (\mathbf{r}) = \mathcal{M} \Phi(\mathbf{r})
   * \f}
   *
   * where the angular flux is a vector containing all the directions, and the
   * moments flux is a vector containing all the moments, i.e.,
   *
   * \f{align*}{
   * & \Psi (\mathbf{r}) \equiv
   * \left[ \begin{array}{c}
   * \Psi_{0} (\mathbf{r}) \\
   * \vdots \\
   * \Psi_{n} (\mathbf{r}) \\
   * \vdots \\
   * \Psi_{N-1} (\mathbf{r}) \\
   * \end{array} \right], \quad
   * & \Phi(\mathbf{r}) \equiv
   * \left[ \begin{array}{c}
   * \Phi_{0}^{0} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{l}^{m} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{L}^{L} (\mathbf{r}) \\
   * \end{array} \right], \quad
   * \f}
   *
   * And the operator \f$ M \f$ is defined as obtaining the projection for a
   * particular direction as follows:
   *
   * \f{align*}{
   * \Psi_n (\mathbf{r}) = \mathcal{M}_{n} \Phi(\mathbf{r})
   * & = \sum_{l = 0}^{L}\sum_{m=-l}^{l}
   * Y_{l}^{m}(\Omega_n) \Phi_{l}^{m}(\mathbf{r}) \\
   * & = \sum_{i = 0}^{M} Y_{i}(\Omega_n) \Phi_{i}(\mathbf{r})
   * & n = 0, \dots, N-1 .
   * \f}
   *
   * In practice we do not want to deal with two indices for the spherical
   * harmonics, so that a unique index is defined by the @ref MomentIndexer
   * class.
   *
   * @note This operator acts independently of the energy group. Thus, when
   * we want to apply this operation to the vector composed by all the energy
   * groups, the operator should be applied independently to each group.
   */
  template <int dim>
  class MomentsToDirections
  {
  public:

    /**
     * @brief Constructor.
     * @param ord
     * @param moments_order
     */
    MomentsToDirections (const std::shared_ptr<QuadratureBase<dim> > & ord,
                         const unsigned int moments_order);

    /**
     * @brief Apply \f$\mathcal{M}\f$
     * @details Given the vector with the angular moments, return the vector
     * with all the directions
     *
     * @f[
     * \Psi (\mathbf{r}) = \mathcal{M} \Phi (\mathbf{r})
     * @f]
     *
     * where each direction has been calculated by applying \f$\mathcal{M}_n\f$
     *
     * @param psi
     * @param phi
     */
    void
    vmult (BlockVector<double> & psi,
           const BlockVector<double> & phi);

    /**
     * @brief It applies \f$\mathcal{M}_n\f$ to the moments.
     * @details It returns the vector \f$ \Psi_n (\mathbf{r}) \f$ calculated by
     *
     * \f{align*}{
     * \Psi_n (\mathbf{r})
     * & = \mathcal{M}_{n} \Phi(\mathbf{r}) \\
     * & = \sum_{i = 0}^{M}
     * Y_{i}(\Omega_n) \Phi_{i}(\mathbf{r})
     * \f}
     *
     * where \f$ i = \text{index}(m,l) \f$ is the moment indexer function to
     * align all indices in a vector to use just one loop.
     * @param psi_n
     * @param phi
     * @param n
     */
    void
    vmult (Vector<double> & psi_n,
           const BlockVector<double> & phi,
           const unsigned int n);

    /**
     * @brief It applies \f$\mathcal{M}_{i \to n}\f$ to the \f$ i\f$-moment.
     * @details It returns the vector \f$ \Psi_n (\mathbf{r}) \f$ calculated by
     *
     * \f{align*}{
     * \Psi_n (\mathbf{r})
     * & = \Psi_n + \mathcal{M}_{i \to n} \Phi_{i}(\mathbf{r}) \\
     * & = \Psi_n + Y_{i}(\Omega_n) \Phi_{i}(\mathbf{r})
     * \f}
     *
     * where \f$ i = \text{index}(m,l) \f$ is the moment indexer function to
     * align all indices in a vector to use just one loop.
     * @note For the zero order moment we perform a direct assignment
     * independent of the particular direction
     * \f{align*}{
     * \Psi_n (\mathbf{r})
     * = \mathcal{M}_{0 \to n} \Phi_{0}(\mathbf{r})
     * = \Phi_{0}(\mathbf{r})
     * \f}
     * @param psi_n
     * @param phi_i
     * @param n
     * @param i
     */
    void
    vmult_add (Vector<double> & psi_n,
               const Vector<double> & phi_i,
               const unsigned int n,
               const unsigned int i);

    /**
     * @brief It applies \f$\mathcal{M}_{0 \to n}\f$ to the \f$ 0\f$-moment.
     * @details It returns the vector \f$ \Psi_n (\mathbf{r}) \f$ calculated by
     *
     * \f{align*}{
     * \Psi_n (\mathbf{r})
     * = \mathcal{M}_{0 \to n} \Phi_{0}(\mathbf{r})
     * = \Phi_{0}(\mathbf{r})
     * \f}
     *
     * where \f$ \Phi_{0} \f$ is the scalar flux, and then the zero
     * order spherical harmonic function is a constant equal to 1, and we
     * have the angular representation independent of the particular direction.
     * @param psi_n
     * @param phi_i
     */
    void
    vmult0 (Vector<double> & psi_n,
            const Vector<double> & phi_i);

    /**
     * @brief Return the memory consumption for an object of this class.
     * @todo Add the memory by the MomentIndexer class.
     */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      return memory_consumption;
    }

  private:

    /** @todo Document me */
    const std::shared_ptr<QuadratureBase<dim> > & mp_ord;

    /** @todo Document me */
    const unsigned int m_n_angles;

    /** @todo Document me */
    const unsigned int m_moments_order;

    /** @todo Document me */
    MomentIndexer<dim> m_moment_indexer;

    /** @todo Document me */
    const unsigned int m_n_moments;
  };

} // end of namespace Forest

#endif /* FOREST_MOMENTS_TO_DIRECTIONS_H */
