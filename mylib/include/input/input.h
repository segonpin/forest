/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input.h
 * @brief  Input class template declarations
 */

#ifndef FOREST_INPUT_H
#define FOREST_INPUT_H

#include "input/input_geom.h"
#include "input/input_mat.h"
#include "input/input_settings.h"

#include <iostream>
#include <string>
#include <map>

/**
 @brief Classes related to data input.
 @defgroup ForestInput ForestInput
 */
namespace Forest
{
  /**
   * @brief Classes related to data input.
   * @details This class acts like a container for the different classes
   * containing the input information
   * @ingroup ForestInput ForestInput
   * @todo add check to the input.
   */
  class Input
  {
  public:

    Input ();

    Input (const std::string & input_file);

    /** @todo document me */
    InputSettings mp_settings;
    /** @todo document me */
    InputGeom mp_geom;
    /** @todo document me */
    InputMat mp_mat;

  public:

    void
    set_input_files (std::string path,
                     std::string filename_no_ext)
    {
      mp_settings.set_files (path, filename_no_ext);
    }

    void
    set_problem (std::string type,
                 bool product_quadrature,
                 std::string quad,
                 std::vector<unsigned int> order)
    {
      mp_settings.set_problem (type, product_quadrature, quad, order);
    }

    void
    set_algebra (bool matrix_free,
                 const bool use_fission_density,
                 const std::string eig_solver,
                 const double eig_tol,
                 const unsigned int eig_max_it,
                 const std::string mg_solver,
                 const double mg_tol,
                 const unsigned int mg_max_it,
                 const std::string inner_solver,
                 const double inner_tol,
                 const unsigned int inner_max_it)
    {
      mp_settings.set_algebra (matrix_free, use_fission_density,
          eig_solver, eig_tol, eig_max_it,
          mg_solver, mg_tol, mg_max_it,
          inner_solver, inner_tol, inner_max_it);
    }

    void
    set_fe_settings (unsigned int degree,
                     unsigned int n_ref)
    {
      mp_settings.set_fe_settings (degree, n_ref);
    }

    void
    set_output (std::vector<bool> & options)
    {
      mp_settings.set_output (options);
    }

    void
    set_geom (unsigned int dim,
              Core core,
              std::map<unsigned int, Lattice> lattices,
              std::map<unsigned int, Pin> pins)
    {
      mp_geom.set(dim, core, lattices, pins);
    }

    void
    set_xs (const unsigned int id,
         const std::vector<double> & st,
         const std::vector<std::vector<double> > & ss,
         const std::vector<double> & nsf,
         const std::vector<double> & chi,
         const std::string name = "Mixture")
    {
      mp_mat.set(id, st, ss, nsf, chi, name);
    };


    void
    check ()
    {
      std::cout << "check geometry" << '\n';
      mp_geom.check();
      std::cout << "check materials" << '\n';
      mp_mat.check();
      std::cout << "check finished" << '\n';
    }

    std::string
    show ()
    {
      return mp_geom.show();
    }
  };

} // end of namespace Forest

#endif /* FOREST_INPUT_H */

