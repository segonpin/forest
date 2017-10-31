#ifndef FOREST_UTILS_TIMER_H
#define FOREST_UTILS_TIMER_H

#include <map>                          // for map
#include <string>                       // for string

namespace Forest
{
  /**
   * @brief StatsTimer
   * @ingroup ForestUtils
   * @note Don't forget to declare copy constructor and assignment
   * operator but don't implement them. You want to make sure they
   * are inaccessible otherwise you may accidentally get copies of
   * your singleton appearing.
   */
  class StatsTimer
  {
  public:

    /** @brief Instance */
    static StatsTimer &
    instance ()
    {
      static StatsTimer _instance;
      return _instance;
    }

    /** add */
    void
    add (const std::string & flag, double value);

    /** add */
    void
    add (const char* flag, double value);

    /** add_final */
    void
    add_final (double value);

    /** print_detailed */
    void
    print_detailed ();

    /** get_total */
    double
    get_total () const;

  private:
    /** @brief Singleton constructor. */
    StatsTimer ()
        : total_time (0),
          final_time (0)
    {
    }

    /** @brief Don't Implement copy constructor for singleton. */
    StatsTimer (StatsTimer const&);

    /** @brief Don't implement assignment operator for singleton. */
    void
    operator= (StatsTimer const&);

    /* Member variables. */
    double total_time;
    double final_time;
    std::map<double, std::string> detailed;
  };

} // end of namespace Forest

#endif
