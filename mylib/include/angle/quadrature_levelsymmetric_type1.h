/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_levelsymmetric_type1.h
 * @brief   QuadratureLevelSymType1 class template declarations
 */

#ifndef FOREST_QUADRATURE_LEVELSYMMETRIC_TYPE1_H
#define FOREST_QUADRATURE_LEVELSYMMETRIC_TYPE1_H

#include "angle/quadrature_base.h"

#include <vector>

//
// Classes used for quadrature. The group is documented in forest_angle.h
//

namespace Forest
{
  /**
   * @brief  QuadratureLevelSymType1
   * @ingroup ForestAngle
   */

  template <int dim>
  class QuadratureLevelSymType1 : public QuadratureBase<dim>
  {
  public:
    QuadratureLevelSymType1 (unsigned int order_);

    static unsigned int
    get_n_angles_per_octant (unsigned int a_order_)
    {
      if (dim == 1)
      {
        return 0;
      }
      else if (dim == 2)
      {
        return ((a_order_ / 2) * (a_order_ / 2 + 1)) / 2;
      }
      else if (dim == 3)
      {
        return ((a_order_ / 2) * (a_order_ / 2 + 1)) / 2;
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
    set_ind_quad (std::vector<double> &cosines,
                  std::vector<double> &tmp_w,
                  unsigned int n);

    void
    get_weight_dist (std::vector<unsigned int> &w_dist) const;

    double
    get_dimension_weight () const
    {
      return ((dim == 1) ? 0.5 : (dim == 2) ? 0.25 : 0.125);
    }

  };
} // end of namespace Forest

#endif
