/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_gausslegendre.h
 * @brief   QuadratureGaussLegendre class template declarations
 */

#ifndef FOREST_GENERATE_GAUSS_CHEBYSHEV_H
#define FOREST_GENERATE_GAUSS_CHEBYSHEV_H

#include "utils/constants.h"

#include <vector>
#include <cmath>

namespace Forest
{

  /**
   * @brief generate_gauss_chebyshev
   * @ingroup ForestAngle
   * Generate the weight of the Gauss-Chebyshev quadrature. Because
   * the quadrature is symmetric, we only calculate half of the abscissas and
   * weights, more specifically the "positive" abscissas
   *
   * @param x
   * @param w
   * @param m
   * @param normalize
   *
   */
  template <int dim>
  inline
  void
  generate_gauss_chebyshev (std::vector<double> &x, std::vector<double> &w,
                            unsigned int m, bool normalize = false)
  {
    // Loop over the desired roots.
    double W = Forest::pi / (double) m;
    for (unsigned int i = 1; i <= m; i++)
    {
      x[i - 1] = std::cos ((i - 0.5) * Forest::pi / m);
    }

    // Generate weights
    for (unsigned int i = 0; i < m; ++i)
    {
      w[i] = W * std::sqrt (1.0 - x[i] * x[i]);
    }

    // Optionally normalize the weights to 2
    if (normalize)
    {
      double w_tot = 0.0;
      for (unsigned int i = 0; i < m; ++i)
      {
        w_tot += w[i];
      }
      for (unsigned int i = 0; i < m; ++i)
      {
        w[i] *= 2.0 / w_tot;
      }
    }
  }

} // end of namespace Forest

#endif
