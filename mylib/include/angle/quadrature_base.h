/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_base.h
 * @brief   QuadratureBase class template declarations
 */

#ifndef FOREST_QUADRATURE_BASE_H
#define FOREST_QUADRATURE_BASE_H

#include <vector>
#include <string>

//
// Classes used for quadrature. The namespace is documented in forest_angle.h
//

/**
 *
 * @todo Have quadratures based on the "octant" and
 * "angle within octant" for the moment, but also general
 * quadratures without this restriction, so I have to develop
 * extra interfaces and a flag for quadratures not defined by octant.
 *
 */

namespace Forest
{
  /**
   * @brief QuadratureBase
   * @ingroup ForestAngle
   */
  template <int dim>
  class QuadratureBase
  {
  public:

    QuadratureBase (const unsigned int n_angles_per_octant_,
                    const std::string & name_);

  protected:

    unsigned int a_order;
    unsigned int p_order;
    unsigned int order;

    unsigned int n_angles_per_octant;
    std::string name;
    unsigned int n_octants;
    unsigned int n_angles;
    unsigned int n_leg_moments;

    std::vector<std::vector<double> > octant_sign;

    std::vector<std::vector<double> > xx;
    std::vector<double> ww;

    std::vector<std::vector<unsigned int> > in2out;

    /** @todo document me */
    void
    set_n_octants ();

    /** @todo document me */
    void
    set_n_angles ();

    /** @todo document me */
    void
    set_octant_sign ();

    /** @todo document me */
    void
    set_all_angles (double angular_weight);

    /** @todo document me */
    void
    set_bc_map ();

    /** @todo document me */
    void
    set_in_from_out ();

    /** @todo document me */
    void
    set_a_order (unsigned int a_order_)
    {
      a_order = a_order_;
    }

    /** @todo document me */
    void
    set_p_order (unsigned int p_order_)
    {
      p_order = p_order_;
    }
    void
    set_order (unsigned int order_)
    {
      order = order_;
    }

    /** @todo document me */
    void
    set_n_leg_moments (unsigned int order_)
    {
      n_leg_moments = order_ / 2;
    }

  public:

    std::vector<std::vector<unsigned int> > out2in;

    /** @todo document me */
    unsigned int
    get_n_angles () const
    {
      return n_angles;
    }

    /** @todo document me */
    unsigned int
    get_n_octants () const
    {
      return n_octants;
    }

    /** @todo document me */
    unsigned int
    get_n_angles_per_octant () const
    {
      return n_angles_per_octant;
    }

    /** @todo document me */
    unsigned int
    get_a_order () const
    {
      return a_order;
    }

    /** @todo document me */
    unsigned int
    get_p_order () const
    {
      return p_order;
    }

    /** @todo document me */
    unsigned int
    get_order () const
    {
      return order;
    }

    /** @todo document me */
    unsigned int
    get_n_leg_moments () const
    {
      return n_leg_moments;
    }

    /**
     *  @brief Return quadrature points \f$ q(:) \f$.
     */
    const std::vector<std::vector<double> > &
    get_q () const
    {
      return xx;
    }

    /**
     *  @brief Return quadrature point \f$ q \f$.
     *  @param a   Angle position
     */
    const std::vector<double> &
    get_q (unsigned int a) const
    {
      return xx[a];
    }

    /**
     *  @brief Return quadrature point \f$ q(d) \f$.
     *  @param a    Angle within octant
     *  @param d    which coordinate
     */
    double
    get_q (unsigned int a,
           unsigned int d) const
    {
      return xx[a][d];
    }

    /**
     *  @brief Return the weights.
     */
    const std::vector<double> &
    get_w () const
    {
      return ww;
    }

    /**
     *  @brief Return the weight.
     *  @param i Angle within octant
     */
    double
    get_w (unsigned int i) const
    {
      return ww[i];
    }

    /**
     * @brief Returns the index of the outgoing direction associated with the
     * incoming direction @p in_angle_ind at the face @p face
     *
     * @param face_
     * @param in_angle_ind_
     * @return out_angle_ind_
     */
    unsigned int
    get_out_ind (unsigned int face_,
                 unsigned int in_angle_ind_) const
    {
      return in2out[face_][in_angle_ind_];
    }

    /** @todo document me */
    double
    get_octant_sign (unsigned int o,
                     unsigned int coordinate) const
    {
      return octant_sign[o][coordinate];
    }

    /**
     *  @brief Return cardinal angle index.
     *  @param o    Octant index
     *  @param a    Angle within octant
     */
    unsigned int
    get_index (unsigned int o,
               unsigned int a) const
    {
      return o * n_angles_per_octant + a;
    }

    /**
     *  @brief Are the indices valid?
     *  @param o    Octant index
     *  @param a    Angle within octant
     */
    bool
    valid_index (unsigned int o,
                 unsigned int a) const
    {
      return (o < n_octants) and (a < n_angles_per_octant);
    }

    // Display and check

    /** @todo document me */
    void
    disp () const;

    /** @todo document me */
    void
    check_ordinates () const;

    /** @todo document me */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      memory_consumption += sizeof(xx) + sizeof(ww) + sizeof(octant_sign)
                            + sizeof(in2out);
      return memory_consumption;
    }

  };
} // end of namespace Forest

#endif
