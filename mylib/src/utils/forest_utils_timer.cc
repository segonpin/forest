/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   utils/forest_utils_timer.cc
 * @brief  Implementation of class Timer
 */

#include "utils/forest_utils_timer.h"
#include "utils/forest_utils_logstream.h"

#include <map>                      // for map
#include <string>                   // for string
#include <vector>                   // for vector
#include <iostream>

//
// Here we have different functions that can be useful for different small
// tasks. The namespace is documented in forest_utils.h
//

namespace Forest
{

  void
  StatsTimer::add_final (const double value)
  {
    final_time = value;
  }

  void
  StatsTimer::add (const std::string & flag, const double value)
  {
    detailed.insert (std::pair<double, std::string> (value, flag));
    total_time += value;
  }

  void
  StatsTimer::add (const char * flag, const double value)
  {
    detailed.insert (
        std::pair<double, std::string> (value, std::string (flag)));
    total_time += value;
  }

  void
  StatsTimer::print_detailed ()
  {
    log_print_line ();
    log_print_title ("Timing statistics (Wall time)");
    log_print_line (' ');
    std::map<double, std::string>::reverse_iterator it;
    for (it = detailed.rbegin (); it != detailed.rend (); ++it)
    {
      log_print_time (it->second, it->first);
    }
    log_print_line (' ');
    log_print_time ("Total time calculated:", total_time);
    if (final_time != 0)
    {
      log_print_line (' ');
      log_print_time ("Final time provided:", final_time);
    }
  }

  double
  StatsTimer::get_total () const
  {
    return total_time;
  }

} // end of namespace Forest
