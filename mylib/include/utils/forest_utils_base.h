#ifndef FOREST_UTILS_BASE_H
#define FOREST_UTILS_BASE_H

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <boost/algorithm/string.hpp> // for is_any_of, split, trim, trim_if

#include <cstdlib>                      // for atoi, atof
#include <iomanip>                      // for operator<<, setw
#include <fstream>
#include <string>                       // for string, allocator, etc
#include <vector>                       // for vector
#include <sstream>
#include <iterator>

namespace Forest
{

  //-------------------------------------------------------------------------
  /* with this we can avoid some boost that generates many warnings */

  /**
   @brief Split string into tokens
   @details Note that this solution does not skip
   empty tokens, so the following will find 4 items,
   one of which is empty:
   @code{.cpp}
   std::vector<std::string> x;
   split("one:two::three", ':', x);
   @endcode
   */
  inline
  void split(const std::string &s, char delim, std::vector<std::string> &elems) {
      std::stringstream ss;
      ss.str(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
        elems.push_back(item);
      }
  }

  /**
   @brief Split string into tokens
   @details Note that this solution does not skip
   empty tokens, so the following will find 4 items,
   one of which is empty:
   @code{.cpp}
   std::vector<std::string> x = split("one:two::three", ':');
   @endcode
   */
  inline
  std::vector<std::string>
  split(const std::string &s, char delim) {
      std::vector<std::string> elems;
      split(s, delim, elems);
      return elems;
  }

  //-------------------------------------------------------------------------
  // File utils

  /**
   @brief Check if the file exist and can open it.
   @param  filename Name of the file we want to open.
   @retval true if the file is open.
   @retval false otherwise.
   */
  inline bool
  file_exists (const std::string & filename)
  {
    std::ifstream input_file (filename.c_str ());
    return input_file.is_open ();
  }

  /*!
   @brief This function look for the path of the folder where the
   file path_plus_file_plus_exts_ is stored.
   @details It can be used as the working directory for the rest of the files.
   If @p path_plus_file_plus_ext_ is only the name of a file, then "."
   is returned as the working directory. For example, using it as follows
   @code{.cpp}
   path = get_path("/this/is/the/path/file.ext1.ext2");
   @endcode
   we obtain that @c path contains "/this/is/the/path/".
   @param path_plus_file_plus_exts_ The full path of the file that we use.
   @retval string The root for the directory.
   @note The character "/" works for Linux, but probably fails for Windows.
   */
  inline std::string
  get_path (const std::string & path_plus_file_plus_exts_)
  {
    std::size_t pos_ = path_plus_file_plus_exts_.find_last_of ("\\/");
    std::string path;
    if (pos_ == path_plus_file_plus_exts_.npos)
      path = std::string ("./");
    else
      path = path_plus_file_plus_exts_.substr (0, pos_) + "/";
    return path;
  }

  /*!
   @brief Get the name of the file in @p path_plus_file_plus_exts_
   stripping the rest of the path.
   @details For example, using it as follows
   @code{.cpp}
   file_plus_ext = get_file_plus_ext("/this/is/the/path/file.ext1.ext2");
   @endcode
   we obtain that @c file_plus_ext contains "file.ext1.ext2".
   @param path_plus_file_plus_ext_ The full path of the file that we use.
   @retval string The name of the input file without the root path.
   */
  inline std::string
  get_file_plus_ext (const std::string & path_plus_file_plus_ext_)
  {
    std::size_t pos_ = path_plus_file_plus_ext_.find_last_of ("\\/");
    std::string file_plus_ext;
    if (pos_ == path_plus_file_plus_ext_.npos)
    {
      file_plus_ext = path_plus_file_plus_ext_;
    }
    else
    {
      file_plus_ext = path_plus_file_plus_ext_.substr (pos_ + 1);
    }
    return file_plus_ext;
  }

  /*!
   @brief This function look for name of the file in the
   path_plus_file_plus_ext_ and strips the last extension (if any).

   @details This function obtains a "problem-related name" that will be used
   later as a root-name to generate the output files obtained with this input.
   For example, using it as follows
   @code{.cpp}
   file_minus_ext = get_file_minus_ext("/this/is/the/path/file.ext1.ext2");
   @endcode
   we obtain that @c file_minus_ext contains
   "/this/is/the/path/file.ext1".

   @param path_plus_file_plus_ext_ The full name of the input file that we use
   @retval string The name of the input without the last extension.
   */
  inline std::string
  get_file_minus_ext (const std::string & path_plus_file_plus_ext_)
  {
    std::size_t pos_ = path_plus_file_plus_ext_.find_last_of (".");
    std::string file_minus_ext;
    if (pos_ == path_plus_file_plus_ext_.npos)
      file_minus_ext = path_plus_file_plus_ext_;
    else
      file_minus_ext = path_plus_file_plus_ext_.substr (0, pos_);
    return file_minus_ext;
  }

  //-------------------------------------------------------------------------
  // String conversion utilities

  /**
   @brief Convert a string to a number for some types
   @param str The string we want to convert to a number
   @retval val The string into a number
   @note The functions were not compiling. To fix it I put them as
   inline functions. Actually I think that they should be inline function
   because the specification of the function is so short and they are
   good candidate for inline functions.
   */
  template <typename Number>
  Number
  string_to_num (std::string str);

  /**
   @brief Convert a string to a @c unsigned @c char
   @param str The string we want to convert to a number
   @retval val The string into a number
   */
  template <>
  inline
  unsigned char
  string_to_num<unsigned char> (std::string str)
  {
    return atoi (str.c_str ());
  }

  /**
   @brief Convert a string to a @c unsigned @c int
   @param str The string we want to convert to a number
   @retval val The string into a number
   */
  template <>
  inline
  unsigned int
  string_to_num<unsigned int> (std::string str)
  {
    return atoi (str.c_str ());
  }

  /**
   @brief Convert a string to a @c int
   @param str The string we want to convert to a number
   @retval val The string into a number
   */
  template <>
  inline
  int
  string_to_num<int> (std::string str)
  {
    return atoi (str.c_str ());
  }

  /**
   @brief Convert a string to a @c double
   @param str The string we want to convert to a number
   @retval val The string into a number
   */
  template <>
  inline
  double
  string_to_num<double> (std::string str)
  {
    return atof (str.c_str ());
  }

  /*!
   @brief Convert a number of different types to a string
   @details The functions were not compiling. To fix it I put them as
   inline functions. Actually I think that they should be inline function
   because the specification of the function is so short and they are
   good candidate for inline functions.
   @param val The number we want to convert to a string
   @retval str The string we want to convert to a number
   */
  template <typename Number>
  inline std::string
  num_to_string (Number val)
  {
    // std::to_string is a C++11 standard
    //return std::to_string(val);

    std::ostringstream ss;
    ss << val;
    return ss.str ();
  }

  //----------------------------------------------------------------------

  /** @todo document me */
  template <typename Number>
  inline
  void
  string_to_vector (const std::string & in,
                    std::vector<Number> & out)
  {
    std::string tmp = in;
    std::vector<std::string> strs;
    boost::trim (tmp);
    //split(in, ',',strs);
    boost::split (strs, tmp, boost::is_any_of (", "), boost::token_compress_on);
    for (unsigned int i = 0; i < strs.size (); ++i)
    {
      out.push_back (string_to_num<Number> (strs[i]));
    }
  }

  /** @todo document me */
  template <typename Number>
  inline
  void
  string_to_vector (const std::string & in,
                    std::vector<std::vector<Number> > & out)
  {
    std::string tmp = in;
    boost::trim_if (tmp, boost::is_any_of ("\n "));
    boost::trim_right_if (tmp, boost::is_any_of (";"));
    std::vector<std::string> strs;
    //split(in, ';',strs);
    boost::split (strs, tmp, boost::is_any_of (";"), boost::token_compress_on);
    out.resize (strs.size ());
    for (unsigned int i = 0; i < strs.size (); ++i)
    {
      string_to_vector (strs[i], out[i]);
    }
  }

  /** @todo document me */
  template <typename Number>
  inline
  void
  vector_to_string (const std::vector<Number> & in,
                    std::string & out,
                    const std::string & indent = std::string (""),
                    const unsigned int n_indent = 0,
                    const unsigned int precision = 7)
  {
    out.clear ();

    // We should check that the vector is not empty
    Assert(in.size () != 0, dealii::ExcMessage("Empty vector!"));

    // We put the numbers inside the string    
    std::stringstream ss;
    //ss << std::fixed;
    ss << std::scientific;
    ss.precision (precision);
    // ss.setf( std::ios::fixed, std::ios::floatfield);

    for (unsigned int i = 0; i < in.size (); ++i)
    {
      ss << std::setw (10) << in[i] << " ";
    }

    //std::copy(in.begin(), in.end(), std::ostream_iterator<Number>(ss, " "));
    for (unsigned int j = 0; j < n_indent; ++j)
      out += indent;

    out += ss.str ();
    out = out.substr (0, out.length () - 1);

  }

  /** @todo document me */
  // We only use the indent for the pretty output
  template <typename Number>
  inline
  void
  vector_to_string (const std::vector<std::vector<Number> > & in,
                    std::string & out,
                    const std::string & indent = std::string (""),
                    const unsigned int n_indent = 0,
                    const unsigned int precision = 7)
  {
    out.clear ();
    //out += "\n";
    std::string bin;
    for (unsigned int i = 0; i < in.size (); ++i)
    {
      vector_to_string (in[i], bin, indent, n_indent,precision);
      //out += bin + ";\n";
      out += bin + ";"; //
    }
    for (unsigned int j = 0; j < n_indent - 1; ++j)
      out += indent;
  }

  //----------------------------------------------------------------------
  // User defined data structures

  /**
   @brief Structure encoding the shape of the reactor
   @details This structure is used to see the shape of the reactor, i.e.,
   \code
   for (unsigned int i = xy_shape.x_begin[], i < xy_shape.x_end[], ++i)
   \endcode
   where the index is for the position along the \c 'y' axis of the \c 'x' line.
   Example of reactor, where \c 'X' means that there is an assembly,
   and \c '0' means that there is not assembly.
   @verbatim
   0 1 2 3 4 - POSITION OVER X AXIS

   0   0 0 X 0 0
   1   0 X X X 0
   2   X X X X X
   3   0 X X X 0
   4   0 0 X 0 0

   x_begin[0] = 2; x_end[0] = 3;
   x_begin[1] = 1; x_end[1] = 4;
   x_begin[2] = 0; x_end[2] = 5;
   x_begin[3] = 1; x_end[3] = 4;
   x_begin[4] = 2; x_end[4] = 3;
   @endverbatim
   */
  typedef struct
  {
    /** index of assembly beginning */
    std::vector<unsigned int> x_begin;
    /** index of assembly end */
    std::vector<unsigned int> x_end;
  } xy_shape;

  //-------------------------------------------------------------------------
  /**
   @brief parsing lines from a file
   */
  template <typename Number>
  void
  parse_multiline_to_vector (std::ifstream &input_file,
                             unsigned int nlines,
                             std::vector<Number> &out,
                             bool deprecate_first = true);

  /**
   @brief parsing lines from a file
   */
  template <typename Number>
  void
  parse_shapeline_to_vector (std::ifstream &input_file,
                             std::vector<unsigned int> &n_nodes,
                             std::vector<xy_shape> &flat_shape,
                             std::vector<Number> &out,
                             bool deprecate_first = true);

  //-------------------------------------------------------------------------

  /**
   @brief Print @p vec to a file @p filename with an @p introduction and @p precision.

   @details It starts to print at the end of the file, adding the vector to
   the previous content of the file.

   @param vec Name of the vector we want to print.
   @param filename Name of the file where we want to print.
   @param introduction The string we want to print before we print the
   vector values. By default it jump two lines.
   @param precision The number of digits we want to print. By default
   it is set equal to 8.
   */
  template <typename Number>
  void
  add_vec_to_file (const std::vector<Number> & vec,
                   const std::string & filename,
                   std::string & introduction = "\n \n",
                   unsigned int precision = 8);

  //-------------------------------------------------------------------------

  /**
   @brief
   @details The notation is following the following convention (for dim == 3):

   @code{.cpp}
   r_normal(0) = {-1,  0,  0}
   r_normal(1) = { 1,  0,  0}
   r_normal(2) = { 0, -1,  0}
   r_normal(3) = { 0,  1,  0}
   r_normal(4) = { 0,  0, -1}
   r_normal(5) = { 0,  0,  1}
   @endcode

   @param f
   @return
   */
  template <int dim>
  inline std::vector<double>
  r_normal (unsigned int f)
  {
    Assert(f < dim * 2, dealii::ExcMessage("Face out of bounds"));
    std::vector<double> r_normal (dim, 0.0);
    if (f % 2 == 0)
      r_normal[(int) f / 2] = -1.;
    else
      r_normal[(int) f / 2] = 1.;
    return r_normal;
  }

} // end of namespace Forest
#endif
