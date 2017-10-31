/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_gausslegendre.h
 * @brief   QuadratureGaussLegendre class template declarations
 */

#ifndef FOREST_QUADRATURE_GAUSS_CHEBYSHEV_H
#define FOREST_QUADRATURE_GAUSS_CHEBYSHEV_H

#include "angle/quadrature_base.h"
#include "utils/constants.h"

#include <vector>

namespace Forest
{
  /**
   * @brief  QuadratureGaussChebyshev
   * @ingroup ForestAngle
   */
  template <int dim>
  class QuadratureGaussChebyshev : public QuadratureBase<dim>
  {
  public:
    QuadratureGaussChebyshev<dim> (unsigned int a_order_,
                                   unsigned int p_order_);

    static unsigned int
    get_n_angles_per_octant (unsigned int a_order_, unsigned int /* p_order_ */)
    {
      if (dim == 1)
      {
        return 0;
      }
      else if (dim == 2 or dim == 3)
      {
        return ((a_order_ + 1) * (a_order_ + 1));
      }
      else
      {
        return 0;
      }
    }

  private:

    unsigned int n_cosines;

    void
    set_number_angles ();

    void
    set_ordinates ();

    void
    set_ind_quad (std::vector<double> &cosines, std::vector<double> &tmp_w,
                  unsigned int n);

    void
    get_weight_dist (std::vector<unsigned int> &w_dist) const;

    double
    get_dimension_weight () const
    {
      return ((dim == 1) ? 0.5 :
              (dim == 2) ? Forest::inv_four_pi : Forest::inv_four_pi);
    }

  };
} // end of namespace Forest

#endif
