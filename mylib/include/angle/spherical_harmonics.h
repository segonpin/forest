/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/spherical_harmonics.h
 * @brief   SphericalHarmonics class template declarations
 */

#ifndef FOREST_SPHERICAL_HARMONICS_H
#define FOREST_SPHERICAL_HARMONICS_H

#include <vector>

namespace Forest
{
  /**
   *  @note Use this class with two arguments: index and direction,
   *  where depending on the dimension the indices will be numbers or vectors
   *  @class SphericalHarmonics
   *  @ingroup ForestAngle
   *  @brief Spherical harmonics generation for anisotropic scattering.
   *
   *  In 2-d and 3-d problems, an expansion of the angular flux in spherical
   *  harmonics is required to treat anisotropic scattering.  In 1-d problems,
   *  the angular flux is expanded in Legendre polynomials, which
   *  can be considered to be a subset of the spherical harmonics.
   *
   *  We follow the presentation of Hebert.
   *
   *  # Spherical Harmonics in 1D (Legendre polynomials).
   *
   *  A function of a directional cosine \f$\xi\f$, \f$ f(\xi) \f$, can
   *  be represented as in \f$ L\f$-th order Legendre expansion:
   *  \f[
   *      f( \xi) \approx \sum^L_{l=0} \frac{2l+1}{2}P_l(\xi)f_l \, ,
   *  \f]
   *  where the Legendre coefficients are defined
   *  \f[
   *      f_l = \int^{1}_{-1} d\xi P_l(\xi) f(\xi) \, .
   *  \f]
   *  In the standard 1-d formulation, in which the azimuthal dependence is
   *  removed by integration, the zeroth Legendre moment angular
   *  flux \f$ \psi(z,\xi) \f$ is
   *  \f[
   *      \phi_0(z) = \int^{1}_{-1} d\xi  \psi(z,\xi) \, ,
   *  \f]
   *  which we recognize as the scalar flux.  The scattering source is
   *  defined in 1D as
   *  \f[
   *      Q(z,\xi) = \int^{2\pi}_{0} d\phi' \int^{1}_{-1}d\xi'
   *                      \Sigma_s(z,\mu_0)\psi(z,\xi') \\
   *  \f]
   * and can likewise be expanded (skipping several important
   * steps! see e.g. the 22106 notes) as
   *  @f[
   *     Q(z,\xi) \approx \sum^{L}_{l=0}
   *     \frac{2l+1}{2}\Sigma_{sl}(z)P_l(\xi)\phi_l(z) \, ,
   *  @f]
   * where the Legendre moments \f$ \Sigma_{sl} \f$ are provided in
   * a cross-section library.  Note, for consistency with the
   * multidimensional coordinates, we place 1D problems along
   * the \f$ z\f$-axis.  For 1D problems using an \f$ L\f$-th order
   * expansion (for scattering), there are \f$ L+1\f$ flux moments
   * per node (and possibly more than one node per mesh).
   *
   *  # Spherical Harmonics in 2D/3D .
   *
   * For 2D and 3D problems, the azimuthal dependence requires use of
   * the spherical harmonics for expansions.
   *
   * The spherical harmonics \f$ Y^m_l(\Omega) = Y^m_l(\xi,\varphi)\f$
   * (in the source, denoted \f$ Y_{lm}\f$ ) are defined
   *  @f[
   *     Y^m_l(\xi,\varphi) = \sqrt{ (2-\delta_{m,0})
   *                                    \frac{(l-|m|)!}{(l+|m|)!} }
   *                             P^{|m|}_l (\xi) \mathcal{T}_m (\varphi ) \ .
   *  @f]
   * where the associated Legendre polynomials are defined in terms of
   * Legendre polynomials as
   *  @f[
   *     P^m_l(\xi) = (1-\xi^2)^{m/2} \frac{d^m}{d\xi^m}
   *                   P_l(\xi) \, , \,\,\,\,\, m \geq 0 \, ,
   *  @f]
   * and the trigonometric functions are
   *  @f[
   *     \mathcal{T}_m (\varphi ) \left\{
   *       \begin{array}{l l}
   *          \cos(m\varphi)      & \quad \text{if $m \geq 0$ }\\
   *          \sin(|m|\varphi)    & \quad \text{otherwise} \, . \\
   *       \end{array} \right.
   *  @f]
   * The associated Legendre polynomials are as given use the so-called Ferrer
   * definition, omitting the typically standard factor of \f$ (-1)^m \f$ .
   * Hebert notes this representation is helpful since
   *  @f[
   *       \mathbf{\Omega} =  \begin{pmatrix}
   *      \mu      \\
   *      \eta     \\
   *      \xi
   *   \end{pmatrix} =  \begin{pmatrix}
   *      \sin(\theta)\cos(\varphi)      \\
   *      \sin(\theta)\sin(\varphi)      \\
   *      \cos(\theta)
   *   \end{pmatrix} = \begin{pmatrix}
   *      Y^{1}_{1}(\xi,\varphi)   \\
   *      Y^{-1}_{1}(\xi,\varphi)  \\
   *      Y^{0}_{1}(\xi,\varphi)
   *     \end{pmatrix} \, .
   *  @f]
   * The spherical harmonics through \f$ L=3\f$ in terms of the directional
   * cosines are
   *  @f[
   *     Y^{0}_{0} = 1 \, ,
   *  @f]
   *  @f[
   *     Y^{-1}_{1}=\eta \, , \quad
   *     Y^{0}_{1}=\xi   \, , \quad
   *     Y^{1}_{1}=\mu
   *  @f]
   *  @f[
   *     Y^{-2}_{2} = \sqrt{3}\mu\eta          \, , \quad
   *     Y^{-1}_{2}=\sqrt{3}\xi\eta            \, , \quad
   *     Y^{0}_{2}=\frac{1}{2}(3\xi^2-1)       \, , \quad
   *     Y^{1}_{2}=\sqrt{3}\xi\mu              \, , \quad
   *     Y^{2}_{2}=\frac{\sqrt{3}}{2}(\mu^2-\eta^2) \, .
   *  @f]
   *  @f[
   *     Y^{-3}_{3} = \sqrt{\frac{5}{8}}\eta(3\mu^2-\eta^2)      \, , \quad
   *     Y^{-2}_{3} = \sqrt{15}\xi\mu\eta                        \, , \quad
   *     Y^{-1}_{3} = \sqrt{\frac{3}{8}}\eta(5\xi^2-1)           \, , \quad
   *     Y^{0}_{3}  = \frac{1}{2}(5\xi^3-3\xi)                   \, , \quad
   *  @f]
   *  @f[
   *     Y^{1}_{3}  = \sqrt{\frac{3}{8}}\mu(5\xi^2-1)            \, , \quad
   *     Y^{2}_{3}  = \sqrt{\frac{15}{4}}\xi(\mu^2-\eta^2)       \, , \quad
   *     Y^{3}_{3}  = \sqrt{\frac{5}{8}}\mu(\mu^2-3\eta^2)       \, .
   *  @f]
   * The trigonometric functions, associated Legendre polynomials, and
   * spherical harmonics satisfy the following orthogonality conditions:
   *  @f[
   *      \int^{\pi}_{-\pi} d\varphi \mathcal{T}_n (\varphi )
   *         \mathcal{T}_{n'}(\varphi ) = \pi(1+\delta_{n,0})\delta_{n,n'} \, ,
   *  @f]
   *  @f[
   *      \int^{1}_{-1} d\mu P^m_l (\mu)
   *        P^{m'}_{l'} = \frac{2(l+m)!}{(2l+1)(l-m)!}\delta_{l,l'} \, ,
   *  @f]
   * and
   *  @f[
   *       \int_{4\pi} d^2\Omega Y^m_l(\mathbf{\Omega})Y^m_{l'}(\mathbf{\Omega})
   *           = \int^{2\pi}_{0} d\phi' \int^{1}_{-1}d\mu'
   *             Y^m_l(\theta,\varphi) Y^{m'}_{l'}(\theta,\varphi)
   *           = \frac{4\pi}{2l+1} \delta_{l,l'}\delta_{m,m'} \, .
   *  @f]
   * Then, angular flux is expanded following
   *  @f[
   *     \psi(\mathbf{r},\mathbf{\Omega}) \approx
   *       \sum^{L}_{l=0} \frac{2l+1}{4\pi}
   *       \sum^{l}_{m=-l} Y^m_l (\mathbf{\Omega}) \phi^m_l(\mathbf{r}) \, ,
   *  @f]
   * where
   *  @f[
   *     \phi^m_l(\mathbf{r}) = \int_{4\pi} d^2\Omega Y^m_l(\mathbf{\Omega})
   *                            \psi(\mathbf{r},\mathbf{\Omega}) \, .
   *  @f]
   * The additional \f$ 2\pi\f$ in the expansion comes from the isotropy in
   * the azimuthal angle, multiplied away in 1D (a 1D angular flux is really
   * in units of 1/rad, not 1/steradian).
   *
   * Again, we note that the zeroth order moment is the angular flux.  For
   * an \f$ L\f$-th order expansion, there are \f$(L+1)^2\f$ moments.
   *
   * @note
   * For 2-D, the moments \f$\phi^0_{l>0}\f$
   * vanish identically.  This is because a 2D problem is defined such
   * that there is no variation in the polar (\f$ z\f$) direction.  The moments
   * \f$\phi^0_{l>0}\f$ represent net changes in the \f$ z\f$ direction, and
   * so should be eliminated for efficiency.  For 2D problems, the number of
   * moments is therefore \f$(L+1)^2-L \f$ . This is handled in
   * @ref MomentsToDirections.
   *
   * It should be noted that the for \f$m=0 \f$, the spherical harmonics as
   * defined reduce to the Legendre polynomials in \f$\xi\f$.  Hence, this
   * class can be used for 1D expansions.
   *
   * Currently, only the spherical harmonics through \f$ L=3 \f$ are
   * implemented, though adding capability for arbitrary orders should be
   * straightforward if external libraries are used.
   *
   * @note
   * Legendre expansions are typically referred to in terms of the "order",
   * which we have denoted via \em L.  In many mathematics texts, the
   * term "degree" is used instead.  Here, both "Legendre order" and
   * "spherical harmonic degree" will refer to the \em l subscript, while
   * the \em m subscript represents the "spherical harmonic order".
   *
   *  Refs
   * - Alain Hebert, <em>Applied Reactor Physics</em>, Presses Internationales
   *   Polytechnique, Montreal, 2009.
   *
   */
  template <int dim>
  class SphericalHarmonics
  {
  public:

    /**
     * @brief Calculate \f$ Y^m_l(\Omega)\f$ given \f$ l,m\f$ and \f$\Omega\f$ .
     * @param lm    Legendre order (or spherical harmonic degree)
     * @param omega Direction for evaluation.
     * @return \f$ Y_l^m(\Omega) \f$
     * @note This is the only public interface for calculating the spherical
     * harmonics on a particular point.
     */
    static double
    Y_lm (const std::vector<int> & lm,
          const std::vector<double> & omega);

    /* Possible test to be added to check orthogonality
     static bool
     test ();
     */

  private:

    /**
     * @brief Calculate \f$ Y^m_l \f$ given \f$\xi\f$ and \f$\varphi\f$ .
     * @note I do not want this function to be called, so it is private.
     * @param l      Legendre order (or spherical harmonic degree)
     * @param m      Spherical harmonic order
     * @param xi     direction cosine w/r to z axis
     * @param varphi azimuthal angle (as defined from x axis)
     * @return       \f$ Y^m_l(\Omega) \f$
     */
    static double
    Y_lm (int l,
          int m,
          double xi,
          double varphi);
    /**
     * @brief Calculate \f$ Y^m_l \f$ given \f$\mathbf{\Omega}\f$.
     * @details This private function could be used to interface with
     * external libraries, leaving the rest of the interfaces
     * and their implementations unchanged.
     * @param l   Legendre order (or spherical harmonic degree)
     * @param m   Spherical harmonic order
     * @param mu  direction cosine w/r to x axis
     * @param eta direction cosine w/r to y axis
     * @param xi  direction cosine w/r to z axis
     * @return    \f$ Y^m_l(\Omega) \f$
     * @todo add proper flags to be used with boost when available
     */
    static double
    get_Y_lm (int l,
              int m,
              double mu,
              double eta,
              double xi);

  };

} // end of namespace Forest

#endif
