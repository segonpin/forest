/**
 * @brief   Generate the weight of the Gauss-Legendre quadrature.
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/generate_gauss_legendre.h
 */

#ifndef FOREST_GENERATE_GAUSS_LEGENDRE_H
#define FOREST_GENERATE_GAUSS_LEGENDRE_H

#include "utils/constants.h"

#include <vector>
#include <cmath>

namespace Forest
{

  /**
   * @details Generate the weight of the Gauss-Legendre quadrature. Because
   * @ingroup ForestAngle
   * the quadrature is symmetric, we only calculate half of the abscissas and
   * weights, more specifically the "positive" abscissas.
   * @param x
   * @param w
   * @param m
   * @param x1
   * @param x2
   */
  template <int dim>
  inline
  void
  generate_gauss_legendre (std::vector<double> &x, std::vector<double> &w,
                           unsigned int m, double x1 = -1.0, double x2 = 1.0)
  {
    // High precision is a good idea for this routine.
    const double eps = 1.0e-14;

    // The roots are symmetric in the interval, so we only have
    // to find half of them.
    const long double xm = 0.5 * (x2 + x1);
    const long double xl = 0.5 * (x2 - x1);

    // Loop over the desired roots.
    unsigned int n = 2 * m;
    for (unsigned int i = 1; i <= m; i++)
    {
      long double z = std::cos (Forest::pi * (i - 0.25) / (n + 0.5));
      long double pp, p1;
      // Starting with the above approximation to the i'th root,
      // we enter the main loop of refinement by Newton's method.
      do
      {
        // Loop up the recurrence relation to get the
        // Legendre polynomial evaluated at z.
        long double p2;
        p1 = 1.0;
        p2 = 0.0;
        for (unsigned int j = 0; j < n; j++)
        {
          long double p3;
          p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
        }
        // p1 is now the desired Legendre polynomial. We next
        // compute pp, its derivative, by a standard relation
        // involving also p2, the polynomial of one lower order.
        pp = n * (z * p1 - p2) / (z * z - 1.0);
        // Newton's method.
        z = z - p1 / pp;
      }
      while (std::abs (p1 / pp) > eps);
      // Scale the root to the desired interval, and put in its
      // symmetric counterpart,  (we invert the order to fit in
      // a natural way with the boundary order for the discrete ordinates).

      x[i - 1] = xm + xl * z;
      // Compute the weight and its symmetric counterpart.
      w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    }

  }

} // end of namespace Forest

#endif
