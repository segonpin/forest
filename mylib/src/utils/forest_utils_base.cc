/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   utils/forest_utils_base.cc
 * @brief  Implementation of utility functions
 */

#include "utils/forest_utils_base.h"

#include <fstream>

namespace Forest
{

  //-------------------------------------------------------------------------

  template <typename Number>
  void
  parse_multiline_to_vector (std::ifstream &input_file, unsigned int nlines,
                             std::vector<Number> &out, bool deprecate_first)
  {
    std::string line;
    std::string bin;
    Number num;

    for (unsigned int j = 0; j < nlines; j++)
    {
      getline (input_file, line);
      std::istringstream iss (line);

      // Reading all the elements for a line (but sometimes the first one).
      bool first_element = true;
      while (iss >> bin)
      {
        if (deprecate_first and first_element)
        {
          first_element = false;
          if ((bin == "#") or (bin == "!") or (bin == "//") or (bin == ""))
          {
            /* If the line is one of previous, we move to the next line. */
            j--;
            continue;
          }
        }
        else
        {
          num = string_to_num<Number> (bin);
          out.push_back (num);
        }
      }
    }
  }

#ifndef DOXYGEN
  // Explicit instantiation
  template void
  parse_multiline_to_vector<unsigned char> (std::ifstream &input_file,
                                            unsigned int nlines,
                                            std::vector<unsigned char> &out,
                                            bool deprecate_first);

  template void
  parse_multiline_to_vector<double> (std::ifstream &input_file,
                                     unsigned int nlines,
                                     std::vector<double> &out,
                                     bool deprecate_first);
#endif

  template <typename Number>
  void
  parse_shapeline_to_vector (std::ifstream &input_file,
                             std::vector<unsigned int> &n_nodes,
                             std::vector<xy_shape> &flat_shape,
                             std::vector<Number> &out, bool deprecate_first)
  {
    std::string line;
    std::string bin;
    xy_shape flat_shape_aux;

    Number num;
    out.clear ();

    for (unsigned int j = 0; j < n_nodes[1]; j++)
    {
      getline (input_file, line);
      std::istringstream iss (line);

      // Reading all the elements for a line (but sometimes the first one).
      unsigned int i = 0;
      bool first_element = true;
      while (iss >> bin)
      {
        if (deprecate_first and first_element)
        {
          first_element = false;
          if ((bin == "#") or (bin == "!") or (bin == "//") or (bin == ""))
          {
            // If the line found is a comment or empty line,
            // we move forward to the next line.
            j--;
            continue;
          }
        }
        else
        {
          num = string_to_num<Number> (bin);
          out.push_back (num);
          ++i;
        }
      }

      unsigned int temp = (n_nodes[0] - i) / 2;
      flat_shape_aux.x_begin.push_back (temp);
      flat_shape_aux.x_end.push_back (n_nodes[0] - temp);
    }
    // Verifying that the dimensions are the good ones
    assert(flat_shape_aux.x_begin.size () == n_nodes[1]);
    assert(flat_shape_aux.x_end.size () == n_nodes[1]);
    flat_shape.push_back (flat_shape_aux);
  }

#ifndef DOXYGEN
  // Explicit instantiation
  template void
  parse_shapeline_to_vector<unsigned char> (std::ifstream &input_file,
                                            std::vector<unsigned int> &n_nodes,
                                            std::vector<xy_shape> &flat_shape,
                                            std::vector<unsigned char> &out,
                                            bool deprecate_first);
#endif

  //-------------------------------------------------------------------------

  template <typename Number>
  void
  add_vec_to_file (const std::vector<Number> & vec,
                   const std::string & filename,
                   std::string & introduction,
                   unsigned int precision)
  {
    std::ofstream out (filename.c_str (), std::ios::app);
    out.precision (precision);
    out << introduction;
    for (unsigned int i = 0; i < vec.size (); i++)
    {
      out << vec[i] << " \n";
    }
    out << "\n";
    out.close ();
  }

//-------------------------------------------------------------------------

} // end of namespace Forest
