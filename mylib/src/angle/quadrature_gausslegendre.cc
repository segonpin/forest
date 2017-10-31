/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_gausslegendre.cc
 * @brief  Implementation of class template QuadratureGaussLegendre
 */

#include "angle/quadrature_gausslegendre.h"

#include "utils/constants.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <cmath>                        // for pow, sqrt, abs, cos

namespace Forest
{
  template <int dim>
  QuadratureGaussLegendre<dim>::QuadratureGaussLegendre (
      const unsigned int p_order_)
      : QuadratureBase<dim> (get_n_angles_per_octant (p_order_),
          "GaussLegendre")
  {
    set_ordinates ();
    this->set_bc_map ();
    this->set_in_from_out ();
  }

  /**
   * @brief
   * Given the lower and upper limits of integration x1 and x2,
   * and given n, this routine returns arrays x[1..n] and w[1..n]
   * of length n, containing the abscissas and weights of the
   * Gauss-Legendre n-point quadrature formula in reverse order
   * (from right to left).
   */
  template <>
  void
  QuadratureGaussLegendre<1>::gaulegxw (std::vector<std::vector<double> > &x,
                                        std::vector<double> &w, unsigned int n,
                                        const double x1, const double x2)
  {
    // High precision is a good idea for this routine.
    const double eps = 1.0e-14;

    // The roots are symmetric in the interval, so we only have
    // to find half of them.
    const long double xm = 0.5 * (x2 + x1);
    const long double xl = 0.5 * (x2 - x1);

    // Loop over the desired roots.
    unsigned int m = (n + 1) / 2;
    for (unsigned int i = 1; i <= m; i++)
    {
      long double z = std::cos (Forest::pi * (i - 0.25) / (n + 0.5));
      // Starting with the above approximation to the i'th root,
      // we enter the main loop of refinement by Newton's method.
      long double pp, p1, p3;
      do
      {
        // Loop up the recurrence relation to get the
        // Legendre polynomial evaluated at z.
        long double p2;
        p1 = 1.0;
        p2 = 0.0;
        for (unsigned int j = 0; j < n; j++)
        {
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
      double x_aux1 = xm - xl * z;
      double x_aux2 = xm + xl * z;

      x[n - i][0] = x_aux1;
      x[i - 1][0] = x_aux2;
      // Compute the weight and its symmetric counterpart.
      double w_aux = 2.0 * xl / ((1.0 - z * z) * pp * pp);
      w[n - i] = w_aux;
      w[i - 1] = w_aux;
    }
  }

  template <int dim>
  void
  QuadratureGaussLegendre<dim>::set_ordinates ()
  {
    Assert(false, dealii::ExcMessage("Not implemented"));
  }

  template <>
  void
  QuadratureGaussLegendre<1>::set_ordinates ()
  {
    gaulegxw (xx, ww, n_angles);
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      ww[i] *= 0.5; // scaling weights
    }
  }

  /// Explicit instantiation of class template QuadratureGaussLegendre<dim>
  template class QuadratureGaussLegendre<1> ;
  template class QuadratureGaussLegendre<2> ;
  template class QuadratureGaussLegendre<3> ;

} // end of namespace Forest
