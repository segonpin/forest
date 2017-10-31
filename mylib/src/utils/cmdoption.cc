/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   utils/cmdoption.cc
 * @brief  Implementation of class CmdOption
 */

#include "utils/cmdoption.h"

#include <algorithm>
#include <cstring>
#include <iostream>

namespace Forest
{

  bool
  CmdOption::exists (char** begin,
                     char** end,
                     const std::string& option)
  {
    return std::find (begin, end, option) != end;
  }

  char*
  CmdOption::get (char ** begin,
                  char ** end,
                  const std::string & option)
  {
    char ** itr = std::find (begin, end, option);
    if (itr != end && ++itr != end)
    {
      return *itr;
    }
    return 0;
  }

} // end of namespace Forest
