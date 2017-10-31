/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/moment_indexer.h
 * @brief   MomentIndexer class declarations
 */

#ifndef FOREST_MOMENT_INDEXER_H
#define FOREST_MOMENT_INDEXER_H

#include <vector>

//
// Classes used for quadrature. The group is documented in forest_angle.h
//
namespace Forest
{

  /**
   *  @class MomentIndexer
   *  @ingroup ForestAngle
   *  @brief Indexes spherical harmonics moments
   *
   *  The one group (or within-group) scattering source in 3-d is defined
   *  @f[
   *      Q(\mathbf{r},\mathbf{\Omega}) =
   *      \sum^L_{l=0} \frac{2l+1}{4\pi}
   *      \sum^l_{m=-l} \Sigma_{sl}(\mathbf{r}) \phi^m_l(\mathbf{r}) .
   *  @f]
   *  and for the representation of the angular flux in Spherical Harmonics
   *  @f[
   *      \Psi(\mathbf{r},\mathbf{\Omega}) =
   *      \sum^L_{l=0} \sum^l_{m=-l}
   *      Y_{l}^{m}(\mathbf{r}) \Phi^m_l(\mathbf{r}) .
   *  @f]
   *
   *  Such sums over \f$l\f$ and \f$m\f$ are commonplace, in generating
   *  sources and also in simply accessing moments sequentially.  This class
   *  builds a vector of \f$(l,m)\f$ pairs that can be used in a single loop
   *  over all flux moments \f$\Phi^m_l\f$,
   *  @f[
   *      \Psi(\mathbf{r},\mathbf{\Omega}) \equiv
   *      \sum^M_{i=0} Y_{i}(\mathbf{r}) \Phi_{i}(\mathbf{r})
   *  @f]
   *  where the sumation is actually done using the re-indexing
   *  \f$ i = \text{ind}(m,l) \f$ and the convention
   *  \f{align}{
   *      & Y_{i}(\mathbf{r}) := Y_{l}^{m}(\mathbf{r}) \\
   *      & \Phi_{i}(\mathbf{r}) := \Phi_{l}^{m}(\mathbf{r}) .
   *  \f}
   *
   *  For 3-d, the moments are ordered as
   *  @f[
   *      [(0,0)] \, , \, \, \, [(1,-1),\, (1,0)\, (1,1)] \, ,
   *      \, \, \, [(2,-2),\,(2,-1)\, \ldots
   *  @f]
   *  For an angular order of \f$L\f$, there are \f$(L+1)^2\f$ moments.
   *
   *  For 2-d, several moments can be eliminated by symmetry.  Only those
   *  moments for which \f$l+m\f$ is even are retained (see Hebert).  Thus,
   *  the moments are ordered as
   *  @f[
   *      [(0,0)] \, , \, \, \, [(1,-1),\, (1,1)] \, ,
   *      \, \, \, [(2,-2),\,(2,0)\, \ldots
   *  @f]
   *  For an angular order of \f$L\f$, there are \f$(L+1)(L+2)/2\f$ moments.
   *
   *  For 1-d, \f$m=0\f$, so that the moments are ordered by \f$l\f$ alone,
   *  from 0 to the given order.  Hence, there are \f$L+1\f$ moments.
   *
   */

  template <int dim>
  class MomentIndexer
  {

  public:

    /**
     *  @brief Constructor.
     *  @param angular_order     Spherical harmonics order of the \em flux.
     */
    MomentIndexer (unsigned int angular_order);

    /**
     *  @brief Return vector of values of m for a given l.
     *  @param  l  Angular degree.
     *  @return m values
     */
    const std::vector<int> &
    get_m_index (unsigned int l) const;

    /**
     * @brief Return angular order.
     * @return
     */
    unsigned int
    get_angular_order () const
    {
      return m_angular_order;
    }

    /// Return number of moments.
    unsigned int
    get_number_of_moments () const
    {
      return m_number_of_moments;
    }

    /**
     *  @brief  Returns l value.
     *  @param  i Cardinal moment index.
     *  @return l value.
     */
    int
    get_l (unsigned int i) const
    {
      return m_l_vec[i];
    }

    /**
     *  @brief Return m value.
     *  @param i Cardinal moment index.
     *  @return m value.
     */
    int
    get_m (unsigned int i) const
    {
      return m_m_vec[i];
    }

    /**
     *  @brief  Returns moment index.
     *  @param  l   l value.
     *  @param  m   m value.
     *  @return Cardinal moment index.
     */
    unsigned int
    get_index (unsigned int l, int m) const
    {
      return m_m_index[l][m];
    }

    /** Print the indices. */
    void
    display () const;

  private:

    /** Legendre order of the \em flux. */
    unsigned int m_angular_order;
    /** Number of moments. */
    unsigned int m_number_of_moments;
    /** Vector of proper \f$m\f$ values for a given \f$l\f$. */
    std::vector<std::vector<int> > m_m_index;
    /** Vector of l values. */
    std::vector<int> m_l_vec;
    /** Vector of m values. */
    std::vector<int> m_m_vec;

    void
    construct ();

  };

} // end of namespace Forest

#endif /* FOREST_MOMENT_INDEXER_H */
