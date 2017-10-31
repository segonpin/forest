/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_chebyshev_legendre.cc
 * @brief  Implementation of class template QuadratureChebyshevLegendre
 */

#include "angle/quadrature_chebyshev_legendre.h"

#include "utils/constants.h"
#include "angle/generate_gauss_legendre.h"
#include "angle/generate_gauss_chebyshev.h"

#include <cmath>                        // for pow, sqrt, abs, cos
#include <cstddef>                      // for size_t
#include <algorithm>                    // for max

namespace Forest
{
  template <int dim>
  QuadratureChebyshevLegendre<dim>::QuadratureChebyshevLegendre (
      const unsigned int a_order_, const unsigned int p_order_)
      : QuadratureBase<dim> (get_n_angles_per_octant (a_order_, p_order_),
            "ChebyshevLegendre"),
        na (a_order_ + 1),
        np (p_order_ + 1),
        phi (na),
        cos_phi (na),
        sin_phi (na),
        azimuth_weight (na),
        theta (np),
        cos_theta (np),
        sin_theta (np),
        polar_weight (np)
  {
    this->set_a_order (a_order_);
    this->set_p_order (p_order_);
    this->set_order (std::max (a_order_, p_order_));
    this->set_n_leg_moments (this->order);

    set_ordinates ();

    this->set_all_angles (get_dimension_weight ());
    this->set_bc_map ();
    this->set_in_from_out ();
  }

  template <int dim>
  void
  QuadratureChebyshevLegendre<dim>::set_ordinates ()
  {
    //------------------------------------------------------------------------//
    // AZIMUTH QUADRATURE
    //------------------------------------------------------------------------//

    // Chebyshev points are just equally spaced phi's in the first quadrant
    {
      // temporary arrays
      //std::vector<double> x (na, 0.0);
      //std::vector<double> w (na, 0.0);
      //generate_gauss_chebyshev (x, w, m, false);

      using Forest::pi;
      for (unsigned int i = 0; i < na; ++i)
      {
        phi[i] = 0.25 * (2 * i + 1) * pi / na;
        cos_phi[i] = std::cos (phi[i]);
        sin_phi[i] = std::sin (phi[i]);
        azimuth_weight[i] = 0.5 * pi / na;
      }
    }
    //------------------------------------------------------------------------//
    // POLAR QUADRATURE
    //------------------------------------------------------------------------//

    {
      // temporary arrays
      std::vector<double> x (np, 0.0);
      std::vector<double> w (np, 0.0);

      // generate parameters
      generate_gauss_legendre<dim> (x, w, np);

      // fill array
      for (unsigned int i = 0; i < np; ++i)
      {
        unsigned int j = np - i - 1;
        cos_theta[j] = x[i];
        sin_theta[j] = std::sqrt (1.0 - x[i] * x[i]);
        polar_weight[j] = w[i];
      }
    }

    //------------------------------------------------------------------------//
    // PRODUCT QUADRATURE
    //------------------------------------------------------------------------//

    double scale = 1.0;
    if (dim == 2)
    {
      scale = 2.0;
    }
    double weight_tot = 0.0;
    size_t n = 0;
    for (unsigned int a = 0; a < na; ++a)
    {
      for (unsigned int p = 0; p < np; ++p, ++n)
      {
        this->xx[n][0] = sin_theta[p] * cos_phi[a];
        this->xx[n][1] = sin_theta[p] * sin_phi[a];
        this->xx[n][2] = cos_theta[p];
        this->ww[n] = scale * polar_weight[p] * azimuth_weight[a];
        weight_tot += this->ww[n];
      } // end polar loop
    } // end azimuth loop
    weight_tot *= this->n_octants;

  }

  /// Explicit instantiation of class template QuadratureChebyshevLegendre<dim>
  template class QuadratureChebyshevLegendre<1> ;
  template class QuadratureChebyshevLegendre<2> ;
  template class QuadratureChebyshevLegendre<3> ;

} // end of namespace Forest
