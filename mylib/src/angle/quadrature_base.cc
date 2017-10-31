/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_base.cc
 * @brief  Implementation of class template QuadratureBase
 */

#include "angle/quadrature_base.h"

#include <cmath>                        // for pow, sqrt, abs, cos
#include <iostream>                     // for cout
#include <string>
#include <vector>
#include <iomanip>

namespace Forest
{
  template <int dim>
  QuadratureBase<dim>::QuadratureBase (const unsigned int n_angles_per_octant_,
                                       const std::string & name_)
      : n_angles_per_octant (n_angles_per_octant_),
        name (name_)
  {
    // setting internal information
    set_n_octants ();
    set_n_angles ();
    set_octant_sign ();

    // resizing the containers for the angles
    const unsigned int max_s_dim = 3;
    xx.resize (n_angles, std::vector<double> (max_s_dim, 0.0));
    ww.resize (n_angles);
  }

  template <int dim>
  void
  QuadratureBase<dim>::set_n_octants ()
  {
    n_octants = std::pow (2, dim);
  }

  template <int dim>
  void
  QuadratureBase<dim>::set_n_angles ()
  {
    n_angles = get_n_angles_per_octant () * get_n_octants ();
  }

  /**
   * @brief Assigning the signs for each octant
   * @details We will run over the eight octants assigning the signs
   * @verbatim
   * axis   x y z
   * o1   = + + +
   * o2   = - + +
   * o3   = - - +
   * o4   = + - +
   * o5   = + + -
   * o6   = - + -
   * o7   = - - -
   * o8   = + - -
   * @endverbatim
   */
  template <int dim>
  void
  QuadratureBase<dim>::set_octant_sign ()
  {
    // allocate the memory for the maximum possible
    const unsigned int max_s_dim = 3;
    const unsigned int max_octants = 8;
    octant_sign.resize (max_octants, std::vector<double> (max_s_dim));
    // first
    octant_sign[0][0] = +1.0;
    octant_sign[0][1] = +1.0;
    octant_sign[0][2] = +1.0;
    // second
    octant_sign[1][0] = -1.0;
    octant_sign[1][1] = +1.0;
    octant_sign[1][2] = +1.0;
    // third
    octant_sign[2][0] = -1.0;
    octant_sign[2][1] = -1.0;
    octant_sign[2][2] = +1.0;
    // fourth
    octant_sign[3][0] = +1.0;
    octant_sign[3][1] = -1.0;
    octant_sign[3][2] = +1.0;
    // fifth
    octant_sign[4][0] = +1.0;
    octant_sign[4][1] = +1.0;
    octant_sign[4][2] = -1.0;
    // sixth
    octant_sign[5][0] = -1.0;
    octant_sign[5][1] = +1.0;
    octant_sign[5][2] = -1.0;
    // seventh
    octant_sign[6][0] = -1.0;
    octant_sign[6][1] = -1.0;
    octant_sign[6][2] = -1.0;
    // eighth
    octant_sign[7][0] = +1.0;
    octant_sign[7][1] = -1.0;
    octant_sign[7][2] = -1.0;
  }

  template <int dim>
  void
  QuadratureBase<dim>::set_all_angles (double dimension_weight)
  {
    for (unsigned int o = 0; o < n_octants; ++o)
    {
      for (unsigned int a = 0; a < n_angles_per_octant; ++a)
      {
        for (unsigned int d = 0; d < dim; ++d)
        {
          xx[get_index (o, a)][d] = octant_sign[o][d] * xx[a][d];
        }
      }
    }

    for (unsigned int a = 0; a < n_angles_per_octant; ++a)
    {
      ww[get_index (0, a)] = ww[a] * dimension_weight;
    }

    for (unsigned int o = 1; o < n_octants; ++o)
    {
      for (unsigned int a = 0; a < n_angles_per_octant; ++a)
      {
        ww[get_index (o, a)] = ww[a];
      }
    }
  }

  /**
   * j_ord = bc_map(f, i_ord), with
   * 
   * \f$ \Omega_{\text{out}} = \Omega_{\text{in}}
   * - 2 n_f (\Omega_{\text{in}} \cdot n_f) \f$
   * 
   * the normals are used in the order specified in 
   * http://www.dealii.org/developer/doxygen/deal.II/structGeometryInfo.html
   * 
   * "normals pointing in -x, x, -y, y, -z, z direction".
   * 
   * @todo Should this function be made independent of dealii? Think about it.
   * 
   */
  template <int dim>
  void
  QuadratureBase<dim>::set_bc_map ()
  {
    const unsigned int faces_per_cell = 2 * dim;
    in2out.resize (faces_per_cell, std::vector<unsigned int> (n_angles, -1));

    // initialize the normal vectors to all the faces
    std::vector<std::vector<double> > normals;
    normals.resize (faces_per_cell, std::vector<double> (dim, 0.));
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      if (dim == 1)
      {
        if (f % 2 == 0)
        {
          normals[f][0] = -1.;
        }
        else
        {
          normals[f][0] = +1.;
        }
      }
      else
      {
        if (f % 2 == 0)
        {
          normals[f][(int) f / dim] = -1.;
        }
        else
        {
          normals[f][(int) f / dim] = +1.;
        }
      }
    }

    const double EPS = 1.e-10;
    std::vector<double> omega_out (dim);
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      for (unsigned int dir_in = 0; dir_in < n_angles; ++dir_in)
      {
        /* We calculate the vector we should look for. */
        double Omega_dot_n = 0.;
        for (unsigned int i = 0; i < dim; ++i)
          Omega_dot_n += xx[dir_in][i] * normals[f][i];

        /* If f is an out-going boundary for i_ord we don't have to calculate
           which one is the incoming flux. */
        if (Omega_dot_n < 0 and std::abs (Omega_dot_n) > EPS)
        {
          /* Calculating the norm of the difference. */
          for (unsigned int i = 0; i < dim; ++i)
          {
            omega_out[i] = xx[dir_in][i] - 2 * normals[f][i] * Omega_dot_n;
          }

          for (unsigned int dir = 0; dir < n_angles; ++dir)
          {
            /* Calculating the norm of the difference */
            double diff_ij = 0;
            for (unsigned int i = 0; i < dim; ++i)
            {
              diff_ij += std::pow (omega_out[i] - xx[dir][i], 2);
            }

            /* If the norm is small enough, we have found the vector
               and stop the j_ord loop. */
            if (diff_ij < EPS)
            {
              in2out[f][dir_in] = dir;
              break;
            }
          }
        }
      }
    }

    bool verbose = false;
    if (verbose)
    {
      std::cout << std::endl << "printing inside set_bc_map" << std::endl;

      for (unsigned int f = 0; f < faces_per_cell; ++f)
      {
        std::cout << "normals[" << f << "] = ";
        for (unsigned int i = 0; i < dim; ++i)
        {
          std::cout << normals[f][i] << " ";
        }
        std::cout << std::endl;
      }

      for (unsigned int i_ord = 0; i_ord < n_angles; ++i_ord)
      {
        std::cout << "xx[" << i_ord << "] = ";
        for (unsigned int i = 0; i < dim; ++i)
        {
          std::cout << xx[i_ord][i] << " ";
        }
        std::cout << std::endl;
      }

      std::cout << "bc_map = \n";
      for (unsigned int f = 0; f < faces_per_cell; ++f)
      {
        std::cout << " face " << f << " = ";
        for (unsigned int i_ord = 0; i_ord < n_angles; ++i_ord)
        {
          std::cout << in2out[f][i_ord] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << std::endl;
    }

  }

  /**
   * j_ord = bc_map(f, i_ord), with
   *
   * \f$ \Omega_{\text{in}} = \Omega_{\text{out}}
   * - 2 n_f (\Omega_{\text{out}} \cdot n_f) \f$
   *
   * the normals are used in the order specified in
   * http://www.dealii.org/developer/doxygen/deal.II/structGeometryInfo.html
   *
   * "normals pointing in -x, x, -y, y, -z, z direction".
   *
   * @todo Should this function be made independent of dealii? Think about it.
   *
   */
  template <int dim>
  void
  QuadratureBase<dim>::set_in_from_out ()
  {
    const unsigned int faces_per_cell = 2 * dim;
    out2in.resize (faces_per_cell, std::vector<unsigned int> (n_angles, -1));

    // initialize the normal vectors to all the faces
    std::vector<std::vector<double> > normals;
    normals.resize (faces_per_cell, std::vector<double> (dim, 0.));
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      if (dim == 1)
      {
        if (f % 2 == 0)
        {
          normals[f][0] = -1.;
        }
        else
        {
          normals[f][0] = +1.;
        }
      }
      else
      {
        if (f % 2 == 0)
        {
          normals[f][(int) f / dim] = -1.;
        }
        else
        {
          normals[f][(int) f / dim] = +1.;
        }
      }
    }

    const double EPS = 1.e-6;
    std::vector<double> omega_in (dim);
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      for (unsigned int dir_out = 0; dir_out < n_angles; ++dir_out)
      {
        // We calculate the vector we should look for.
        double Omega_dot_n = 0.;
        for (unsigned int i = 0; i < dim; ++i)
        {
          Omega_dot_n += xx[dir_out][i] * normals[f][i];
        }

        // if Omega_dot_n > 0 out_dir is outgoing flux at face f
        if (Omega_dot_n > 0 and std::abs (Omega_dot_n) > EPS)
        {
          for (unsigned int i = 0; i < dim; ++i)
          {
            omega_in[i] = xx[dir_out][i] - 2 * normals[f][i] * Omega_dot_n;
          }

          for (unsigned int dir = 0; dir < n_angles; ++dir)
          {
            // calculating the norm of the difference
            double diff_ij = 0;
            for (unsigned int i = 0; i < dim; ++i)
            {
              diff_ij += std::pow (omega_in[i] - xx[dir][i], 2);
            }

            // if the norm is small enough, we have found the vector
            // and stop the in_dir loop
            if (diff_ij < EPS)
            {
              out2in[f][dir_out] = dir;
              break;
            }
          }
        }
      }
    }
  }

  /// Displaying the ordinates and the weights in one dimension
  template <int dim>
  void
  QuadratureBase<dim>::disp () const
  {
    std::cout << "Quadrature <" << name << "> being used:" << std::endl;
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      std::cout << "    angle " << std::setw (3) << i << ": " << "    ww[i] = "
                << std::setw (8) << ww[i] << ", " << "    xx[i] = (";

      for (unsigned int d = 0; d < dim - 1; ++d)
      {
        std::cout << std::setw (9) << xx[i][d] << " , ";
      }
      std::cout << std::setw (9) << xx[i][dim - 1] << "); ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  template <>
  void
  QuadratureBase<1>::check_ordinates () const
  {
    // double tol = 1.e-8;
    // checking if the sum of all weights is equal to 2.0 up to a tolerance of tol
    double sum = 0.0;
    for (unsigned int i = 0; i < this->get_n_angles (); ++i)
    {
      sum += ww[i];
    }
    //Assert_msg(std::abs(sum-2.0) < tol,"Sum of all weights is different from 2.0");
  }

  // Explicit generation of class template QuadratureBase<dim>
  template class QuadratureBase<1> ;
  template class QuadratureBase<2> ;
  template class QuadratureBase<3> ;

} // end of namespace Forest
