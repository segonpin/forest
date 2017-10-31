/**
 *  @author Sebastian Gonzalez-Pintor
 *  @date   04/09/2014
 *  @file   utils/cmdoption.h
 *  @brief  CmdOption class declarations.
 */

#ifndef FOREST_CMDOPTION_H
#define FOREST_CMDOPTION_H

#include <string> // for char_traits, string

namespace Forest
{

  /**
   @brief We use this class to parse the command line arguments.
   @ingroup ForestUtils
   @details  The usage is as follows:
   @code{.cpp}
   int main (int argc, char **argv)
   {
     parse_args.exists (argv, argv + argc, "-f")
     Forest::CmdOption parse_args;
     assert(parse_args.exists (argv, argv + argc, "-f"));
     std::string input_file = parse_args.get (argv, argv + argc, "-f");
     return 0;
   }
   @endcode

   @note Actually, the main point is specifying the input file, and
   getting the path. Information about the function can be found here:
   <a href="http://stackoverflow.com/questions/865668/parse-command-line-arguments">
   stackoverflow</a>

   @note Also <a href="http://www.boost.org/doc/libs/release/libs/program_options/">
   boost </a> can be very useful, especially when parsing parameter files.
   */
  class CmdOption
  {
  public:
    /**
     @brief Empty constructor
     */
    CmdOption ()
    {
    }

    /**
     @brief Check if @c option exists in the command line arguments.

     @param begin
     @param end
     @param option
     @return @c true if the option exists and @c false otherwise.
     */
    bool
    exists (char** begin,
            char** end,
            const std::string & option);

    /**
     @brief Get the @c option from the command line, and return it.

     @param begin
     @param end
     @param option
     @return The argument associated to the @p option specified
     */
    char*
    get (char ** begin,
         char ** end,
         const std::string & option);
  };

} // end of namespace Forest

#endif
