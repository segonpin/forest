#ifndef FOREST_UTILS_LOGSTREAM_H
#define FOREST_UTILS_LOGSTREAM_H

#include "deal.II/base/logstream.h"     // for LogStream, deallog
#include "utils/cmdoption.h"
#include "utils/forest_utils_timer.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>                       // for string

//
// Here we have different functions that can be useful for different small
// tasks. The namespace is documented in forest_utils.h
//

namespace Forest
{

  using namespace dealii;

  /** @brief Line length used by default in all the printing */
  static const unsigned int line_length = 80;

  /** @brief Indent length used by default in all the printing */
  static const unsigned int indent_length = 2;

  /** @brief Print time statistics message */
  inline
  void
  log_print_time (const std::string &msg,
                  double value)
  {
    const std::string prefix = deallog.get_prefix ();

    const std::string label = "TIME  ";
    deallog.pop ();
    deallog.push (prefix + label);

    const std::string units = " sec.";
    const unsigned int time_length = 20;
    const unsigned int pref_length = prefix.size () + 2 + label.size ();

    deallog << std::setw (line_length - time_length - pref_length) << std::left
            << std::setfill ('.') << msg;
    deallog << std::setw (time_length - units.size ()) << std::right
            << std::setfill ('.') << std::scientific << value << units
            << std::flush << std::endl;

    deallog.pop ();
    deallog.push (prefix.substr (0, prefix.size () - 1));
  }

  /** @brief Print memory consumption message */
  inline
  void
  log_print_mem (const std::string &msg,
                 double value)
  {
    const std::string prefix = deallog.get_prefix ();

    const std::string label = "MEMORY";
    deallog.pop ();
    deallog.push (prefix + label);
    const std::string units = "  MB.";
    const double Byte_to_MB = 1. / 1.048576e+6;

    const std::string mem_msg = " memory consumption:";

    const unsigned int mem_length = 20;
    const unsigned int pref_length = prefix.size () + 2 + label.size ();

    deallog << std::setw (line_length - mem_length - pref_length) << std::left
            << std::setfill ('.') << msg + mem_msg;
    deallog << std::setw (mem_length - units.size ()) << std::right
            << std::setfill ('.') << std::scientific << value * Byte_to_MB
            << units << std::flush << std::endl;

    deallog.pop ();
    deallog.push (prefix.substr (0, prefix.size () - 1));
  }

  /**
   * @brief Print statistics if required and return convergence status.
   * @details The expected output is as follows:
   * @verbatim
   FOREST:TABLE::Iter.   0:   keff = 1.15459   spw =   21.0   err = 1.54e-01
   FOREST:TABLE::Iter.   1:   keff = 1.16049   spw =   18.0   err = 5.10e-03
   FOREST:TABLE::Iter.   2:   keff = 1.16184   spw =   15.0   err = 1.16e-03
   FOREST:TABLE::Iter.   3:   keff = 1.16215   spw =   14.0   err = 2.60e-04
   ...
   @endverbatim
   * but it can be tuned by using some parameters.
   */
  inline
  bool
  log_table_keff_wgit (const unsigned int iter,
                       const double keff,
                       const double spg,
                       const double err,
                       const double tol,
                       const unsigned int print_freq = 1,
                       const unsigned int indent_times = 1)
  {
    /* Are we converged? */
    const bool converged = err <= tol;
    /* Do we want to say anything? */
    const bool verbose = iter % print_freq == 0;
    /* If we are converged, or we want to print, we continue. */
    if (converged or verbose)
    {
      /* Define some constants. */
      const unsigned int indent_length = 2;
      const unsigned int space_length = 3;
      const std::string prefix = deallog.get_prefix ();
      const std::string label = "TABLE";
      char sep = ' ';
      std::string indent (indent_length * indent_times, sep);
      std::string space (space_length, sep);
      /* Start the deallog. */
      deallog.pop ();
      deallog.push (prefix + label);
      /* Print the information. */
      deallog << indent;
      deallog << "Iter. " << std::setfill (sep) << std::setw (3) << iter << ":"
              << space;
      deallog << "keff = " << std::setfill (sep) << std::fixed
              << std::setprecision (5) << keff << space;
      deallog << "spg = " << std::setfill (sep) << std::fixed
              << std::setprecision (1) << std::setw (5) << spg << space;
      deallog << "err = " << std::setfill (sep) << std::scientific
              << std::setprecision (2) << err << std::flush << std::endl;
      /* Back to original deallog. */
      deallog.pop ();
      deallog.push (prefix.substr (0, prefix.size () - 1));
    }
    /* return the convergence state. */
    return converged;
  }

  /**
   * @brief Print statistics if required and return convergence status.
   * @details The expected output is as follows:
   * @verbatim
   FOREST:TABLE::  Iter.   0:   keff = 1.15408   err = 1.54e-01
   FOREST:TABLE::  Iter.   1:   keff = 1.16051   err = 5.56e-03
   FOREST:TABLE::  Iter.   2:   keff = 1.16185   err = 1.15e-03
   FOREST:TABLE::  Iter.   3:   keff = 1.16215   err = 2.55e-04
   ...
   @endverbatim
   * but it can be tuned by using some parameters.
   */
  inline
  bool
  log_table_keff (const unsigned int iter,
                  const double keff,
                  const double err,
                  const double tol,
                  const unsigned int print_freq = 1,
                  const unsigned int indent_times = 1)
  {
    /* Are we converged? */
    const bool converged = err <= tol;
    /* Do we want to say anything? */
    const bool verbose = iter % print_freq == 0;
    /* If we are converged, or we want to print, we continue. */
    if (converged or verbose)
    {
      /* Define some constants. */
      const unsigned int indent_length = 2;
      const unsigned int space_length = 3;
      const std::string prefix = deallog.get_prefix ();
      const std::string label = "TABLE";
      char sep = ' ';
      std::string indent (indent_length * indent_times, sep);
      std::string space (space_length, sep);
      /* Start the deallog. */
      deallog.pop ();
      deallog.push (prefix + label);
      /* Print the information. */
      deallog << indent;
      deallog << "Iter. " << std::setfill (sep) << std::setw (3) << iter << ":"
              << space;
      deallog << "keff = " << std::setfill (sep) << std::fixed
              << std::setprecision (5) << keff << space;
      deallog << "err = " << std::setfill (sep) << std::scientific
              << std::setprecision (2) << err << std::flush << std::endl;
      /* Back to original deallog. */
      deallog.pop ();
      deallog.push (prefix.substr (0, prefix.size () - 1));
    }
    /* return the convergence state. */
    return converged;
  }

  /**
   * @brief Print @p msg indenting @p indent_times
   * @param msg
   * @param indent_times
   */
  inline
  void
  log_print_text (const std::string & msg,
                  unsigned int indent_times = 0)
  {
    std::string indent (indent_length * indent_times, ' ');
    const unsigned int pref_length = deallog.get_prefix ().size () + 1;
    deallog
        << indent
        << std::setw (line_length - pref_length - indent_length * indent_times)
        << std::left << std::setfill (' ') << msg << std::flush << std::endl;
  }

  /**
   * @brief Print a line by repeating the character in @p sep
   * @param sep = '-'
   */
  inline
  void
  log_print_line (char sep = '-')
  {
    const unsigned int pref_length = deallog.get_prefix ().size () + 1;
    deallog << std::setw (line_length - pref_length) << std::right
            << std::setfill (sep) << sep << std::flush << std::endl;
  }

  /**
   * @brief Print title centered @p end_msg filling empty space with @p sep
   * @param end_msg
   * @param sep = ' '
   */
  inline
  void
  log_print_title (const std::string & end_msg,
                   char sep = ' ')
  {
    const unsigned int pref_length = deallog.get_prefix ().size () + 1;
    const unsigned int left_length = (line_length - pref_length
                                      - end_msg.size ())
                                     / 2;
    const unsigned int right_length = line_length - pref_length
                                      - end_msg.size () - left_length;
    deallog << std::string (left_length, sep) << end_msg
            << std::string (right_length, sep) << std::flush << std::endl;
  }

  inline std::string
  start_forest (int argc,
                char **argv,
                const unsigned int output_level = 1)
  {
    // Parsing the arguments from the command line
    Forest::CmdOption parse_args;
    Assert(parse_args.exists (argv, argv + argc, "-f"),
        dealii::ExcMessage("Input not specified."));
    std::string input_file = parse_args.get (argv, argv + argc, "-f");
    static std::ofstream ostream_file (input_file + ".log");

    // Preparing the stream for the console output
    deallog.attach (ostream_file);
    deallog.pop ();
    deallog.push ("FOREST");
    deallog.test_mode (false);
    deallog.depth_console (output_level);
    deallog.depth_file (output_level);

    /* Print starting message*/
    log_print_line ('*');
    log_print_title ("FOREST:");
    log_print_title ("(F)ramework (O)riented to (RE)actor (S)imula(T)ions");
    log_print_line ('*');

    /* We return the input file for calling the Input class. */
    return input_file;
  }

  inline
  void
  finish_forest ()
  {
    Forest::StatsTimer::instance ().print_detailed ();
    log_print_line ('*');
    log_print_title ("Calculation finished", '.');
    log_print_line ('*');
  }

} // end of namespace Forest

#endif
