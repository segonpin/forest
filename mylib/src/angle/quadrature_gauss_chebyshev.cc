/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_gausslegendre.cc
 * @brief  Implementation of class template QuadratureGaussLegendre
 */

#include "angle/quadrature_gauss_chebyshev.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <cmath>                        // for pow, sqrt, abs, cos
#include <iostream>                     // for cout
#include <algorithm>
#include <cstdlib>

namespace Forest
{
  template <int dim>
  QuadratureGaussChebyshev<dim>::QuadratureGaussChebyshev (
      const unsigned int a_order_, const unsigned int p_order_)
      : QuadratureBase<dim> (get_n_angles_per_octant (a_order_, p_order_),
          "GaussChebyshev")
  {
    set_ordinates ();

    this->set_a_order (a_order_);
    this->set_p_order (p_order_);
    this->set_order (std::max (a_order_, p_order_));
    this->set_n_leg_moments (this->order);
    this->set_all_angles (get_dimension_weight ());
    this->set_bc_map ();
    this->set_in_from_out ();
  }

  /**
   * @brief Hardcoding several values for the discrete ordinates generation.
   */
  template <int dim>
  void
  QuadratureGaussChebyshev<dim>::set_ind_quad (std::vector<double> &tmp_q,
                                               std::vector<double> &tmp_w,
                                               unsigned int order_)
  {
    double mu;
    switch (order_)
    {
      case 2:
        mu = std::sqrt (1.0 / 3.0);
        tmp_w[0] = 1.0;
        break;

      case 4:
        mu = 0.3500212;
        tmp_w[0] = 1.0 / 3.0;
        break;

      case 6:
        mu = 0.2666355;
        tmp_w[0] = 0.1761263;
        tmp_w[1] = 0.1572071;
        break;

      case 8:
        mu = std::sqrt (1.0 / 21.0);
        tmp_w[0] = 0.1209877;
        tmp_w[1] = 0.0907407;
        tmp_w[2] = 0.0925926;
        break;

      case 12:
        mu = 0.1672126;
        tmp_w[0] = 0.0707626;
        tmp_w[1] = 0.0558811;
        tmp_w[2] = 0.0373377;
        tmp_w[3] = 0.0502819;
        tmp_w[4] = 0.0258513;
        break;

      case 16:
        mu = 0.1389568;
        tmp_w[0] = 0.0489872;
        tmp_w[1] = 0.0413296;
        tmp_w[2] = 0.0212326;
        tmp_w[3] = 0.0256207;
        tmp_w[4] = 0.0360486;
        tmp_w[5] = 0.0144586;
        tmp_w[6] = 0.0344958;
        tmp_w[7] = 0.0085179;
        break;

      default:
        // change this for an error flag at exit the function
        std::cout
            << "Error:\n"
            << "  Quadrature levelsymmetric type 1 is not implemented for SN"
            << std::endl << "  different from 2, 4, 6, 8, 12 and 16."
            << std::endl
            << "  Thus 'abort()' is called to terminate the program execution."
            << std::endl << std::endl;
        abort ();
        break;
    }

    // Number of different points over an axis

    tmp_q[0] = mu;
    double mu_sq = mu * mu;
    for (unsigned int i = 1; i < n_cosines; ++i)
    {
      tmp_q[i] = std::sqrt (mu_sq + i * (1 - 3 * mu_sq) / (n_cosines - 1));
    }
  }

  /**
   *
   * S_2  weight distribution
   *                  {1};
   *
   * S_4  weight distribution
   *                  {1,
   *                  1,1};
   *
   * S_6  weight distribution
   *                  {1,
   *                  2,2,
   *                 1,2,1};
   *
   * S_8  weight distribution
   *                  {1,
   *                  2,2,
   *                 2,3,2,
   *                1,2,2,1};
   *
   * S_10 weight distribution
   *                  {1,
   *                  2,2,
   *                 3,4,3,
   *                2,4,4,2,
   *               1,2,3,2,1};
   *
   * S_12 weight distribution
   *                  {1,
   *                  2,2,
   *                 3,4,3,
   *                3,5,5,3,
   *               2,4,5,4,2,
   *              1,2,3,3,2,1};
   *
   * S_14 weight distribution
   *                  {1,
   *                  2,2,
   *                 3,5,3,
   *                4,6,6,4,
   *               3,6,7,6,3,
   *              2,5,6,6,5,2,
   *             1,2,3,4,3,2,1};
   *
   * S_16 weight distribution
   *                  {1,
   *                  2,2,
   *                 3,5,3,
   *                4,6,6,4,
   *               4,7,8,7,4,
   *              3,6,8,8,6,3,
   *             2,5,6,7,6,5,2,
   *            1,2,3,4,4,3,2,1};
   *
   * S_18 weight distribution
   *                  {1,
   *                 2,  2,
   *               3,  6,  3,
   *             4,  7,  7,  4,
   *           5,  8,  9,  8,  5,
   *         4,  8, 10, 10,  8,  4,
   *       3,  7,  9, 10,  9,  7,  3,
   *     2,  6,  7,  8,  8,  7,  6,  2,
   *   1,  2,  3,  4,  5,  4,  3,  2,  1};
   *
   * S_20 weight distribution
   *                  {1,
   *                 2,  2,
   *               3,  6,  3,
   *             4,  7,  7,  4,
   *           5,  8, 10,  8,  5,
   *         5,  9, 11, 11,  9,  5,
   *       4,  8, 11, 12, 11,  8,  4,
   *     3,  7, 10, 11, 11, 10,  7,  3,
   *   2,  6,  7,  8,  9,  8,  7,  6,  2,
   * 1,  2,  3,  4,  5,  5,  4,  3,  2,  1};
   *
   */
  template <int dim>
  void
  QuadratureGaussChebyshev<dim>::get_weight_dist (
      std::vector<unsigned int> &w_dist) const
  {

    switch (this->order)
    {
      case 2:
        w_dist =
        { 1};
        break;

        case 4:
        w_dist =
        { 1,
          1, 1};
        break;

        case 6:
        w_dist =
        {
          1,
          2, 2,
          1, 2, 1};
        break;

        case 8:
        w_dist =
        {
          1,
          2, 2,
          2, 3, 2,
          1, 2, 2, 1};
        break;

        case 10:
        w_dist =
        {
          1,
          2, 2,
          3, 4, 3,
          2, 4, 4, 2,
          1, 2, 3, 2, 1};
        break;

        case 12:
        w_dist =
        {
          1,
          2, 2,
          3, 4, 3,
          3, 5, 5, 3,
          2, 4, 5, 4, 2,
          1, 2, 3, 3, 2, 1};
        break;

        case 14:
        w_dist =
        {
          1,
          2, 2,
          3, 5, 3,
          4, 6, 6, 4,
          3, 6, 7, 6, 3,
          2, 5, 6, 6, 5, 2,
          1, 2, 3, 4, 3, 2, 1};
        break;

        case 16:
        w_dist =
        {
          1,
          2, 2,
          3, 5, 3,
          4, 6, 6, 4,
          4, 7, 8, 7, 4,
          3, 6, 8, 8, 6, 3,
          2, 5, 6, 7, 6, 5, 2,
          1, 2, 3, 4, 4, 3, 2, 1};
        break;

        case 18:
        w_dist =
        {
          1,
          2, 2,
          3, 6, 3,
          4, 7, 7, 4,
          5, 8, 9, 8, 5,
          4, 8, 10, 10, 8, 4,
          3, 7, 9, 10, 9, 7, 3,
          2, 6, 7, 8, 8, 7, 6, 2,
          1, 2, 3, 4, 5, 4, 3, 2, 1};
        break;

        case 20:
        w_dist =
        {
          1,
          2, 2,
          3, 6, 3,
          4, 7, 7, 4,
          5, 8, 10, 8, 5,
          5, 9, 11, 11, 9, 5,
          4, 8, 11, 12, 11, 8, 4,
          3, 7, 10, 11, 11, 10, 7, 3,
          2, 6, 7, 8, 9, 8, 7, 6, 2,
          1, 2, 3, 4, 5, 5, 4, 3, 2, 1};
        break;
      }
    }

  template <int dim>
  void
  QuadratureGaussChebyshev<dim>::set_ordinates ()
  {
    std::vector<double> x_tmp (n_cosines);
    std::vector<double> w_tmp (n_cosines);
    std::vector<unsigned int> w_dist (this->n_angles_per_octant);

    set_ind_quad (x_tmp, w_tmp, this->order);
    get_weight_dist (w_dist);

    unsigned int k = 0;
    for (unsigned int i = 0; i < n_cosines; ++i)
    {
      for (unsigned int j = 0; j < i + 1; ++j, ++k)
      {
        this->xx[k][0] = x_tmp[j];
        this->xx[k][1] = x_tmp[i - j];
        this->xx[k][2] = std::sqrt (
            1.0 - this->xx[k][0] * this->xx[k][0]
            - this->xx[k][1] * this->xx[k][1]);
        this->ww[k] = w_tmp[w_dist[k] - 1];
      }
    }

    x_tmp.clear ();
    w_dist.clear ();
    w_tmp.clear ();

  }

  template <>
  void
  QuadratureGaussChebyshev<1>::set_ordinates ()
  {
    std::cout << "exception is thrown" << std::endl;
    Assert(false, dealii::ExcMessage("Not implemented."));
  }

  /// Explicit instantiation of class template QuadratureGaussChebyshev<dim>
  template class QuadratureGaussChebyshev<1> ;
  template class QuadratureGaussChebyshev<2> ;
  template class QuadratureGaussChebyshev<3> ;

} // end of namespace Forest
