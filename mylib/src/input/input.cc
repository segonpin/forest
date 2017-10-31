/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input.cc
 * @brief  Implementation of class template Input
 */

#include "input/input.h"

namespace Forest
{

  Input::Input ()
      : mp_settings (),
        mp_geom (),
        mp_mat ()
  {

  }

  /**
   * @brief Constructor.
   * @details The first thing is to take the full_path_input and get the relative path
   * (to the executable folder) of the file. Then we store the name of the file
   * without extension and use it for all the input files.
   * 
   */
  Input::Input (const std::string & input_file)
      : mp_settings (input_file),
        mp_geom (mp_settings.get_file_geom ()),
        mp_mat (mp_settings.get_file_mat ())
  {

  }

} // end namespace Forest

