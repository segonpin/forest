/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_gausslegendre.h
 * @brief   QuadratureGaussLegendre class template declarations
 */

#ifndef FOREST_QUADRATURE_CHEBYSHEV_LEGENDRE_H
#define FOREST_QUADRATURE_CHEBYSHEV_LEGENDRE_H

#include "angle/quadrature_base.h"
#include "utils/constants.h"

#include <vector>

namespace Forest
{
  /**
   * @brief  QuadratureChebyshevLegendre
   * @ingroup ForestAngle
   */
  template <int dim>
  class QuadratureChebyshevLegendre : public QuadratureBase<dim>
  {
  public:
    QuadratureChebyshevLegendre<dim> (unsigned int a_order_,
                                      unsigned int p_order_);

    static unsigned int
    get_n_angles_per_octant (unsigned int a_order_, unsigned int p_order_)
    {
      if (dim == 1)
      {
        return 0;
      }
      else if (dim == 2 or dim == 3)
      {
        return ((a_order_ + 1) * (p_order_ + 1));
      }
      else
      {
        return 0;
      }
    }

  private:

    unsigned int na;
    unsigned int np;

    std::vector<double> phi;
    std::vector<double> cos_phi;
    std::vector<double> sin_phi;
    std::vector<double> azimuth_weight;

    std::vector<double> theta;
    std::vector<double> cos_theta;
    std::vector<double> sin_theta;
    std::vector<double> polar_weight;

    void
    set_ordinates ();

    double
    get_dimension_weight () const
    {
      return (
          (dim == 1) ? 0.5 :
          (dim == 2) ? Forest::inv_four_pi : Forest::inv_four_pi * 0.5);
    }

  };
} // end of namespace Forest

#endif
