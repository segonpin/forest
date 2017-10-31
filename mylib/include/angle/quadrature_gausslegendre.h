/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_gausslegendre.h
 * @brief   QuadratureGaussLegendre class template declarations
 */

#ifndef FOREST_QUADRATURE_GAUSSLEGENDRE_H
#define FOREST_QUADRATURE_GAUSSLEGENDRE_H

#include "angle/quadrature_base.h"

#include <vector>

//
// Classes used for quadrature. The group is documented in forest_angle.h
//

namespace Forest
{
  /**
   * @brief  QuadratureGaussLegendre
   * @ingroup ForestAngle
   */
  template <int dim>
  class QuadratureGaussLegendre : public QuadratureBase<dim>
  {
  public:
    QuadratureGaussLegendre (unsigned int p_order_);

    /**
     * We return the number of angles per octant for the
     * present quadrature. We make this function static
     * so it can be accessed from outside the class,
     * in order to call the base class with the number
     * of angles per octant
     *
     * @param p_order_ polar order
     * @return number of angles per octant
     */
    static unsigned int
    get_n_angles_per_octant (unsigned int p_order_)
    {
      if (dim == 1)
      {
        return p_order_;
      }
      else
      {
        return 0;
      }
    }

  private:

    /** @todo document me */
    void
    set_ordinates ();

    /** @todo document me */
    void
    gaulegxw (std::vector<std::vector<double> > &x,
              std::vector<double> &w,
              unsigned int n,
              double x1 = -1.0,
              double x2 = 1.0);

  };
} // end of namespace Forest

#endif
