/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_levelsymmetric_type2.cc
 * @brief  Implementation of class template QuadratureLevelSymType2
 */

#include "angle/quadrature_levelsymmetric_type2.h"

#include <cmath>                        // for pow, sqrt, abs, cos
#include <cstdlib>                      // for abort
#include <iostream>                     // for cout

// I can not use my dbc functions until I get rid of the dealii includes
//#include "utils/design_by_contract.h"

namespace Forest
{
  template <int dim>
  QuadratureLevelSymType2<dim>::QuadratureLevelSymType2 (
      const unsigned int order_)
      : QuadratureBase<dim> (get_n_angles_per_octant (order_), "LevelSymType2"),
        n_cosines (order_ / 2)
  {
    this->set_a_order (order_);
    this->set_p_order (order_);
    this->set_order (order_);
    this->set_n_leg_moments (order_);

    set_ordinates ();

    this->set_all_angles (get_dimension_weight ());
    this->set_bc_map ();
    this->set_in_from_out ();
  }

  /**
   * @brief Hardcoding several values for the discrete ordinates generation.
   */
  template <int dim>
  void
  QuadratureLevelSymType2<dim>::set_ind_quad (std::vector<double> &tmp_q,
                                              std::vector<double> &tmp_w,
                                              unsigned int order_)
  {
    switch (order_)
    {
      case 2:
      {
        tmp_q[0] = 0.577350269189625764509149;
        tmp_w[0] = 1.000000000000000000000000;
        break;
      }

      case 4:
      {
        tmp_q[0] = 0.350021174581540677777041;
        tmp_q[1] = 0.868890300722201205229788;
        tmp_w[0] = 1.0 / 3.0;
        break;
      }

      case 6:
      {
        tmp_q[0] = 0.266635401516704720331535;
        tmp_q[1] = 0.681507726536546927403750;
        tmp_q[2] = 0.926180935517489107558380;
        tmp_w[0] = 0.176126130863383433783565;
        tmp_w[1] = 0.157207202469949899549768;
        break;
      }

      case 8:
      {
        tmp_q[0] = 0.218217890235992381266097;
        tmp_q[1] = 0.577350269189625764509149;
        tmp_q[2] = 0.786795792469443145800830;
        tmp_q[3] = 0.951189731211341853132399;
        tmp_w[0] = 0.120987654320987654320988;
        tmp_w[1] = 0.0907407407407407407407407;
        tmp_w[2] = 0.0925925925925925925925926;
        break;
      }

      case 10:
      {
        tmp_q[0] = 0.189321326478010476671494;
        tmp_q[1] = 0.508881755582618974382711;
        tmp_q[2] = 0.694318887594384317279217;
        tmp_q[3] = 0.839759962236684758403029;
        tmp_q[4] = 0.963490981110468484701598;
        tmp_w[0] = 0.0893031479843567214704325;
        tmp_w[1] = 0.0725291517123655242296233;
        tmp_w[2] = 0.0450437674364086390490892;
        tmp_w[3] = 0.0539281144878369243545650;
        break;
      }

      case 12:
      {
        tmp_q[0] = 0.167212652822713264084504;
        tmp_q[1] = 0.459547634642594690016761;
        tmp_q[2] = 0.628019096642130901034766;
        tmp_q[3] = 0.760021014833664062877138;
        tmp_q[4] = 0.872270543025721502340662;
        tmp_q[5] = 0.971637719251358378302376;
        tmp_w[0] = 0.0707625899700910439766549;
        tmp_w[1] = 0.0558811015648888075828962;
        tmp_w[2] = 0.0373376737588285824652402;
        tmp_w[3] = 0.0502819010600571181385765;
        tmp_w[4] = 0.0258512916557503911218290;
        break;
      }

      case 14:
      {
        tmp_q[0] = 0.151985861461031912404799;
        tmp_q[1] = 0.422156982304796966896263;
        tmp_q[2] = 0.577350269189625764509149;
        tmp_q[3] = 0.698892086775901338963210;
        tmp_q[4] = 0.802226255231412057244328;
        tmp_q[5] = 0.893691098874356784901111;
        tmp_q[6] = 0.976627152925770351762946;
        tmp_w[0] = 0.0579970408969969964063611;
        tmp_w[1] = 0.0489007976368104874582568;
        tmp_w[2] = 0.0227935342411872473257345;
        tmp_w[3] = 0.0394132005950078294492985;
        tmp_w[4] = 0.0380990861440121712365891;
        tmp_w[5] = 0.0258394076418900119611012;
        tmp_w[6] = 0.00826957997262252825269908;
        break;
      }

      case 16:
      {
        tmp_q[0] = 0.138956875067780344591732;
        tmp_q[1] = 0.392289261444811712294197;
        tmp_q[2] = 0.537096561300879079878296;
        tmp_q[3] = 0.650426450628771770509703;
        tmp_q[4] = 0.746750573614681064580018;
        tmp_q[5] = 0.831996556910044145168291;
        tmp_q[6] = 0.909285500943725291652116;
        tmp_q[7] = 0.980500879011739882135849;
        tmp_w[0] = 0.0489872391580385335008367;
        tmp_w[1] = 0.0413295978698440232405505;
        tmp_w[2] = 0.0203032007393652080748070;
        tmp_w[3] = 0.0265500757813498446015484;
        tmp_w[4] = 0.0379074407956004002099321;
        tmp_w[5] = 0.0135295047786756344371600;
        tmp_w[6] = 0.0326369372026850701318409;
        tmp_w[7] = 0.0103769578385399087825920;
        break;
      }

      case 18:
      {
        tmp_q[0] = 0.129344504545924818514086;
        tmp_q[1] = 0.368043816053393605686086;
        tmp_q[2] = 0.504165151725164054411848;
        tmp_q[3] = 0.610662549934881101060239;
        tmp_q[4] = 0.701166884252161909657019;
        tmp_q[5] = 0.781256199495913171286914;
        tmp_q[6] = 0.853866206691488372341858;
        tmp_q[7] = 0.920768021061018932899055;
        tmp_q[8] = 0.983127661236087115272518;
        tmp_w[0] = 0.0422646448843821748535825;
        tmp_w[1] = 0.0376127473827281471532380;
        tmp_w[2] = 0.0122691351637405931037187;
        tmp_w[3] = 0.0324188352558815048715646;
        tmp_w[4] = 0.00664438614619073823264082;
        tmp_w[5] = 0.0312093838436551370068864;
        tmp_w[6] = 0.0160127252691940275641645;
        tmp_w[7] = 0.0200484595308572875885066;
        tmp_w[8] = 0.000111409402059638628382279;
        tmp_w[9] = 0.0163797038522425240494567;
        break;
      }

      case 20:
      {
        tmp_q[0] = 0.120603343036693597409418;
        tmp_q[1] = 0.347574292315847257336779;
        tmp_q[2] = 0.476519266143665680817278;
        tmp_q[3] = 0.577350269189625764509149;
        tmp_q[4] = 0.663020403653288019308783;
        tmp_q[5] = 0.738822561910371432904974;
        tmp_q[6] = 0.807540401661143067193530;
        tmp_q[7] = 0.870852583760463975580977;
        tmp_q[8] = 0.929863938955324566667817;
        tmp_q[9] = 0.985347485558646574628509;
        tmp_w[0] = 0.0370210490604481342320295;
        tmp_w[1] = 0.0332842165376314841003910;
        tmp_w[2] = 0.0111738965965092519614021;
        tmp_w[3] = 0.0245177476959359285418987;
        tmp_w[4] = 0.0135924329650041789567081;
        tmp_w[5] = 0.0318029065936585971501960;
        tmp_w[6] = 0.00685492401402507781062634;
        tmp_w[7] = 0.0308105481755299327227893;
        tmp_w[8] = -0.000139484716502602877593527;
        tmp_w[9] = 0.00544675187330776223879437;
        tmp_w[10] = 0.00474564692642379971238396;
        tmp_w[11] = 0.0277298541009064049325246;
        break;
      }

      default:
      {
        // change this for an error flag at exit the function
        std::cout
            << "Error:\n" << "  Quadrature levelsymmetric type 2(" << order_
            << ") is not implemented for SN" << std::endl
            << "  different from 2, 4, 6, 8, 10, 12, 14, 16, 18 and 20."
            << std::endl
            << "  Thus 'abort()' is called to terminate the program execution."
            << std::endl << std::endl;
        abort ();
      }
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
  QuadratureLevelSymType2<dim>::get_weight_dist (
      std::vector<unsigned int> &w_dist) const
  {

    switch (this->order)
    {
      case 2:
        w_dist =
        {
          1};
        break;

        case 4:
        w_dist =
        {
          1,
          1,1};
        break;

        case 6:
        w_dist =
        {
          1,
          2,2,
          1,2,1};
        break;

        case 8:
        w_dist =
        {
          1,
          2,2,
          2,3,2,
          1,2,2,1};
        break;

        case 10:
        w_dist =
        {
          1,
          2,2,
          3,4,3,
          2,4,4,2,
          1,2,3,2,1};
        break;

        case 12:
        w_dist =
        {
          1,
          2,2,
          3,4,3,
          3,5,5,3,
          2,4,5,4,2,
          1,2,3,3,2,1};
        break;

        case 14:
        w_dist =
        {
          1,
          2,2,
          3,5,3,
          4,6,6,4,
          3,6,7,6,3,
          2,5,6,6,5,2,
          1,2,3,4,3,2,1};
        break;

        case 16:
        w_dist =
        {
          1,
          2,2,
          3,5,3,
          4,6,6,4,
          4,7,8,7,4,
          3,6,8,8,6,3,
          2,5,6,7,6,5,2,
          1,2,3,4,4,3,2,1};
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
  QuadratureLevelSymType2<dim>::set_ordinates ()
  {
    std::vector<double> x_tmp (n_cosines);
    std::vector<double> w_tmp (this->n_angles_per_octant);
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

    w_dist.clear ();
    x_tmp.clear ();
    w_tmp.clear ();
  }

  template <>
  void
  QuadratureLevelSymType2<1>::set_ordinates ()
  {
    std::cout << "exception is thrown" << std::endl;
    //Assert_msg(false, "ExcNotImplemented");
  }

  // Explicit instantiation of class template QuadratureLevelSymType2<dim>
  template class QuadratureLevelSymType2<1> ;
  template class QuadratureLevelSymType2<2> ;
  template class QuadratureLevelSymType2<3> ;

} // end of namespace Forest
