#ifndef FOREST_UTILS_MEMORY_H
#define FOREST_UTILS_MEMORY_H

#include <map>                          // for map
#include <string>                       // for string

namespace Forest
{
  /** @brief StatsMemory.
   *  @ingroup ForestUtils
   */
  class StatsMemory
  {
  public:

    /** instance */
    static StatsMemory &
    instance ()
    {
      static StatsMemory _instance;
      return _instance;
    }

    /** add */
    void
    add (const std::string & flag, double value);

    /** add */
    void
    add (const char* flag, double value);

    /** print_detailed */
    void
    print_detailed ();

    /** get_total */
    double
    get_total () const;

  private:
    /** @brief Singleton constructor. */
    StatsMemory ()
        : total_memory (0)
    {
    }

    /** @brief Don't Implement copy constructor for singleton. */
    StatsMemory (StatsMemory const&);

    /** @brief Don't implement assignment operator for singleton. */
    void
    operator= (StatsMemory const&);

    /* Member variables. */
    double total_memory;
    std::multimap<double, std::string> detailed;
  };

} // end of namespace Forest

#endif
