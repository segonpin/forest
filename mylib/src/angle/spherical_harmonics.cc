/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/spherical_harmonics.cc
 * @brief  Implementation of class SphericalHarmonics
 */

#include "angle/spherical_harmonics.h"

#include "utils/constants.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

//#ifdef FOREST_ENABLE_BOOST
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
//#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace Forest
{

  /** @todo use the macros to enable or disable boost. */
  template <int dim>
  double
  SphericalHarmonics<dim>::Y_lm (const std::vector<int> & /* lm */,
                                 const std::vector<double> & /* omega */)
  {
    Assert(false, dealii::ExcMessage("Dimension dim not implemented."));
    return 0.;
  }

  template <>
  double
  SphericalHarmonics<1>::Y_lm (const std::vector<int> & lm,
                               const std::vector<double> & omega)
  {
    /* Check the input. */
    Assert(lm.size() == 2, dealii::ExcMessage("Wrong index size."));
    Assert(omega.size() == 3, dealii::ExcMessage("Wrong Omega dimension."));
    /* Indices for spherical harmonics. */
    const int l = lm[0];
    const int m = 0;
    /* Direction. */
    const double mu = 0.;
    const double eta = 0.;
    const double xi = omega[0];
    /* Evaluation the Spherical harmonic in the required direction. */
    return get_Y_lm (l, m, mu, eta, xi);
  }

  template <>
  double
  SphericalHarmonics<2>::Y_lm (const std::vector<int> & lm,
                               const std::vector<double> & omega)
  {
    /* Check the input. */
    Assert(lm.size() == 2, dealii::ExcMessage("Wrong index size."));
    Assert(omega.size() == 3, dealii::ExcMessage("Wrong Omega dimension."));
    /* Indices for spherical harmonics. */
    const int l = lm[0];
    const int m = lm[1];
    /* Direction. */
    const double mu = omega[0];
    const double eta = omega[1];
    const double xi = std::sqrt(1.0 - mu*mu - eta*eta);
    /* Evaluation the Spherical harmonic in the required direction. */
    return get_Y_lm (l, m, mu, eta, xi);
  }

  template <>
  double
  SphericalHarmonics<3>::Y_lm (const std::vector<int> & lm,
                               const std::vector<double> & omega)
  {
    /* Check the input. */
    Assert(lm.size() == 2, dealii::ExcMessage("Wrong index size."));
    Assert(omega.size() == 3, dealii::ExcMessage("Wrong Omega dimension."));
    /* Indices for spherical harmonics. */
    const int l = lm[0];
    const int m = lm[1];
    /* Direction. */
    const double mu = omega[0];
    const double eta = omega[1];
    const double xi = omega[2];
    /* Evaluation the Spherical harmonic in the required direction. */
    return get_Y_lm (l, m, mu, eta, xi);
  }


  /* Possible test to add to check that they are orthogonal
  template <int dim>
  bool
  test ()
  {
    // variables for the angular coarsening
    std::string q_type = (dim == 1 ? "GaussLegendre" : "LevelSymType2");
    unsigned int q_n = 3;
    QuadraturePtr<dim> ord = QuadratureFactory<dim>::build (q_n, q_type);
    std::vector<std::vector<double> > q = ord->get_q ();
    const unsigned int n_angles = ord-> get_n_angles();
    for (unsigned int n = 0; n < n_angles; ++n)
      q[n][2] = std::sqrt (1 - (q[n][0] * q[n][0] + q[n][1] * q[n][1]));
    std::vector<double> w = ord->get_w ();

    MomentIndexer<dim> m_ind (ord-> get_n_order());

    const unsigned int n_moments = m_ind.get_number_of_moments ();

    //Test orthogonality for the spherical harmonics
    for (unsigned int i = 0; i < n_moments; ++i)
    {
      std::cout << "(l,m) = ( " << std::setw (4) << m_ind.get_l (i) << " , "
      << std::setw (4) << m_ind.get_m (i) << " ) = ";
      for (unsigned int ind2 = 0; ind2 < n_moments; ++ind2)
      {
        double val = 0;
        for (unsigned int n = 0; n < n_angles; ++n)
        {
          const std::vector<int> lm =
          { m_ind.get_l (i), m_ind.get_m (i) };
          val += w[n] * SphericalHarmonics<dim>::Y_lm (lm, q[n])
                 * SphericalHarmonics<dim>::Y_lm (lm, q[n])
                 * (2 * m_ind.get_l (i) + 1);
        }
        std::cout << std::setprecision (6) << std::setw (15) << val << " , ";
      }
      std::cout << std::endl;
    }
  }*/

  // Calculate the spherical harmonic of degree l, order m
  // given cos(polar) and azimuthal angle
  template <int dim>
  double
  SphericalHarmonics<dim>::Y_lm (const int l,
                                 const int m,
                                 const double xi,
                                 const double varphi)
  {
    Assert(dim == 3, dealii::ExcMessage("Wrong dimension."));
    //Require(xi >= -1.0);
    //Require(xi <= 1.0);
    //Require(varphi >= 0.0);
    //Require(varphi <= Forest::two_pi);
    double sintheta = std::sqrt (1.0 - xi * xi);
    double mu = sintheta * std::cos (varphi);
    double eta = sintheta * std::sin (varphi);
    return get_Y_lm (l, m, mu, eta, xi);
  }


  template <int dim>
  double
  SphericalHarmonics<dim>::get_Y_lm (const int l,
                                     const int m,
                                     const double mu,
                                     const double eta,
                                     const double xi)
  {
    //Require(l >= 0);
    //Require(m <= l);
    //Require(m >= -l);
    //Require(mu >= -1.0);
    //Require(mu <= 1.0);
    //Require(eta >= -1.0);
    //Require(eta <= 1.0);
    //Require(xi >= -1.0);
    //Require(xi <= 1.0);

    /* If the point is not unit norm, it should be pure z coordinate. */
    double unity = mu * mu + eta * eta + xi * xi;
    if (not (std::fabs (1.0 - unity) < 1.0e-5))
      Assert((std::fabs(mu) < 1.e-6) and (std::fabs(eta) < 1.e-6),
          dealii::ExcMessage("If 3d vector is not unit norm, should be 1d."));

    // if (not (std::fabs (1.0 - unity) < 1.0e-5))
    //   Require( (std::fabs(mu) < 1.e-6) && (std::fabs(eta) < 1.e-6) );

    if (true)
    {
      //double phi = std::acos (mu / std::sqrt (1.0 - xi * xi));
      double phi = std::atan2 (eta, mu);
      double theta = std::acos (xi);
      // Normalization so that phi_00 = phi, phi_1m = current components, etc.
      // This leads to a match with definitions given by Hebert.
      double del = (m == 0) ? 1.0 : 2.0;
      double norm = std::pow (-1.0, m)
          * std::sqrt (del * Forest::four_pi / (2. * l + 1.));
      if (m >= 0)
        return norm * boost::math::spherical_harmonic_r (l, m, theta, phi);
      else
        return norm * boost::math::spherical_harmonic_i (l, -m, theta, phi);
    }

    if (l == 0)
      return 1.0;
    else if (l == 1)
    {
      if (m == -1)
        return eta;
      else if (m == 0)
        return xi;
      else if (m == 1)
        return mu;
    }
    else if (l == 2)
    {
      if (m == -2)
        return 1.732050807568877 * mu * eta;
      else if (m == -1)
        return 1.732050807568877 * xi * eta;
      else if (m == 0)
        return 1.5 * xi * xi - 0.5;
      else if (m == 1)
        return 1.732050807568877 * xi * mu;
      else if (m == 2)
        return 0.8660254037844385 * (mu * mu - eta * eta);
    }
    else if (l == 3)
    {
      if (m == -3)
        return 0.7905694150420948 * eta * (3 * mu * mu - eta * eta);
      else if (m == -2)
        return 3.872983346207417 * xi * mu * eta;
      else if (m == -1)
        return 0.6123724356957945 * eta * (5.0 * xi * xi - 1.0);
      else if (m == 0)
        return 2.5 * xi * xi * xi - 1.5 * xi;
      else if (m == 1)
        return 0.6123724356957945 * mu * (5.0 * xi * xi - 1.0);
      else if (m == 2)
        return 1.936491673103708 * xi * (mu * mu - eta * eta);
      else if (m == 3)
        return 0.7905694150420948 * mu * (mu * mu - 3.0 * eta * eta);
    }
    else
    {
/*#ifdef FOREST_ENABLE_BOOST
      //double phi = std::acos (mu / std::sqrt (1.0 - xi * xi));
      double phi = std::atan2 (eta, mu);
      double theta = std::acos (xi);
      // Normalization so that phi_00 = phi, phi_1m = current components, etc.
      // This leads to a match with definitions given by Hebert.
      double del = (m == 0) ? 1.0 : 2.0;
      double norm = std::pow (-1.0, m) *
          std::sqrt (del * Forest::four_pi / (2.*l+1.));
      if (m >= 0)
      return norm * boost::math::spherical_harmonic_r (l, m, theta, phi);
      else
      return norm * boost::math::spherical_harmonic_i (l, -m, theta, phi);
#else
      // Degree not implemented.
      //Assert_msg(false, 'Maximum Legendre order is 3. For higher orders, enable Boost.');
      return 0;
#endif*/
      Assert(false, dealii::ExcMessage("Index not implemented without BOOST."));
    }
    //std::cout << "FOREST_ENABLE_BOOST = " << FOREST_ENABLE_BOOST << std::endl;
    return 0;
  }

 template class SphericalHarmonics<1> ;
 template class SphericalHarmonics<2> ;
 template class SphericalHarmonics<3> ;

} // end of namespace Forest
