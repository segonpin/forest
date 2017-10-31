/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/directions_to_moments.h
 * @brief   DirectionsToMoments class template declarations
 */

#ifndef FOREST_DIRECTIONS_TO_MOMENTS_H
#define FOREST_DIRECTIONS_TO_MOMENTS_H

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
   * @brief Directions to Moments class
   * @ingroup ForestAngle
   * @details Given a discrete ordinates representation of the angular
   * flux, we project it into a moment expansion of the angular flux in
   * spherical harmonics as follows
   *
   * \f{align*}{
   * \Phi (\mathbf{r}) = \mathcal{D} \Psi(\mathbf{r})
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
   * & \Phi(\mathbf{r})
   * \equiv \left[ \begin{array}{c}
   * \Phi_{0} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{i} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{M} (\mathbf{r}) \\
   * \end{array} \right]
   * \equiv \left[ \begin{array}{c}
   * \Phi_{0}^{0} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{l}^{m} (\mathbf{r}) \\
   * \vdots \\
   * \Phi_{L}^{L} (\mathbf{r}) \\
   * \end{array} \right], \quad
   * \f}
   *
   * @note In practice we do not want to deal with two indices for the spherical
   * harmonics, so that a unique index is defined by the @ref MomentIndexer
   * by mapping \f$ i:=\text{ind}(l,m) \f$ .
   *
   * And the operator \f$ D \f$ is defined as obtaining the projection over a
   * particular moment as follows:
   *
   * \f{align*}{
   * \Phi_{i}(\mathbf{r}) = \mathcal{D}_{i} \Psi (\mathbf{r})
   * := \mathcal{D}_{l}^{m} \Psi (\mathbf{r})
   * & = \int_{4\pi} Y_{l}^{m}(\Omega) \Psi(\mathbf{r},\Omega)
   * \ \text{d} \Omega  \\
   * & \approx \sum_{n=0}^{N-1}
   * \mu_{n} Y_{l}^{m}(\Omega_{n}) \Psi_{n}(\mathbf{r}) \quad
   * l = 1, \dots, L; \quad  -l \geq m \leq l  \\
   * & := \sum_{n=0}^{N-1}
   * \mu_{n} Y_{i}(\Omega_{n}) \Psi_{n}(\mathbf{r}) \quad
   * i = 0, \dots, M .
   * \f}
   *
   * where \f$ \{\Omega_n, \mu_n \}_{n=0}^{N-1} \f$ is the set of quadrature
   * points and weights used for the discrete ordinates representation.
   *
   */
  template <int dim>
  class DirectionsToMoments
  {
  public:

    /**
     * @brief Constructor.
     * @param ord
     * @param moments_order
     */
    DirectionsToMoments (std::shared_ptr<QuadratureBase<dim> > & ord,
                         unsigned int moments_order);

    /**
     * @brief Applies \f$ \mathcal{D} \f$ .
     * @details It applies the  full operator
     *
     * \f{align*}{
     * \Phi (\mathbf{r}) = \mathcal{D} \Psi(\mathbf{r})
     * \f}
     *
     * where each moment is calculated applying \f$\mathcal{D}_{i}\f$ .
     * @param phi
     * @param psi
     */
    void
    vmult (BlockVector<double> & phi,
           const BlockVector<double> & psi);

    /**
     * @brief Applies \f$ \mathcal{D}_{i} \f$ .
     * @details Apply the operator \f$ \mathcal{D}_{i} \f$ defined as follows:
     *
     * \f{align*}{
     * \Phi_{i}(\mathbf{r}) = \mathcal{D}_{i} \Psi (\mathbf{r})
     * := \mathcal{D}_{l}^{m} \Psi (\mathbf{r})
     * & = \int_{4\pi} Y_{l}^{m}(\Omega) \Psi(\mathbf{r},\Omega)
     * \ \text{d} \Omega  \\
     * & \approx \sum_{n=0}^{N-1}
     * \mu_{n} Y_{l}^{m}(\Omega_{n}) \Psi_{n}(\mathbf{r}) \quad
     * l = 1, \dots, L; \quad  -l \geq m \leq l  \\
     * & =: \sum_{n=0}^{N-1}
     * \mu_{n} Y_{i}(\Omega_{n}) \Psi_{n}(\mathbf{r}) \quad
     * i = 0, \dots, M .
     * \f}
     * @param phi_i
     * @param psi
     * @param i
     */
    void
    vmult (Vector<double> & phi_i,
           const BlockVector<double> & psi,
           const unsigned int i);

    /**
     * @brief Add one direction contribution to all moments .
     * @details Add the operator \f$ \mathcal{D}_{n\to i} \f$ for a particular
     * direction \f$ n \f$ to all moments \f$ i = 0, \ldots, M\f$ .
     * @param phi
     * @param psi_n
     * @param n
     */
    void
    vmult_add (BlockVector<double> & phi,
               const Vector<double> & psi_n,
               const unsigned int n);


    /**
     * @brief Perform \f$\Phi_{i}=\Phi_{i}+\mathcal{D}_{i \gets n}\Psi_{n}\f$ .
     * @details Add the operator \f$ \mathcal{D}_{i \gets n} \f$ to the vector
     * of the angular moments, where the addition for each moments is defined by
     *
     * \f{align*}{
     * \Phi_{i}(\mathbf{r})
     * & = \Phi_{i}(\mathbf{r})+\mathcal{D}_{i \gets n} \Psi_{n}(\mathbf{r}) \\
     * & = \Phi_{i}(\mathbf{r})+\mu_{n}Y_{i}(\Omega_{n})\Psi_{n}(\mathbf{r})
     * \f}
     * @param phi_i
     * @param psi_n
     * @param i
     * @param n
     */
    void
    vmult_add (Vector<double> & phi_i,
               const Vector<double> & psi_n,
               const unsigned int i,
               const unsigned int n);

    /**
     * @brief Return the memory consumption of the object.
     * @todo Add MomentIndexer.
     */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      return memory_consumption;
    }

  private:

    /** @todo Document me */
    std::shared_ptr<QuadratureBase<dim> > & mp_ord;

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

#endif /* FOREST_DIRECTIONS_TO_MOMENTS_H */
