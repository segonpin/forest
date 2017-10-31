/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   utils/forest_utils_memory.cc
 * @brief  Implementation of class StatsMemory
 */

#include "utils/forest_utils_memory.h"
#include "utils/forest_utils_logstream.h"

#include <map>                       // for map
#include <utility>                   // for pair
#include <string>                    // for string

//
// Here we have different functions that can be useful for different small
// tasks. The namespace is documented in forest_utils.h
//

namespace Forest
{

  void
  StatsMemory::add (const std::string & flag, const double value)
  {
    detailed.insert (std::pair<double, std::string> (value, flag));
    total_memory += value;
  }

  void
  StatsMemory::add (const char * flag, const double value)
  {
    detailed.insert (
        std::pair<double, std::string> (value, std::string (flag)));
    total_memory += value;
  }

  void
  StatsMemory::print_detailed ()
  {
    log_print_line ();
    log_print_title ("Memory statistics");
    log_print_line (' ');
    std::map<double, std::string>::reverse_iterator it;
    for (it = detailed.rbegin (); it != detailed.rend (); ++it)
    {
      log_print_mem (it->second, it->first);
    }
    log_print_line (' ');
    log_print_mem ("Total memory:", total_memory);
  }

  double
  StatsMemory::get_total () const
  {
    return total_memory;
  }

} // end of namespace Forest
