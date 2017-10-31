#ifndef REAL_ASSEMBLY_H
#define REAL_ASSEMBLY_H

#include "input/input_geom.h"    // for Pin (ptr only), Lattice (ptr only)

#include "deal.II/base/point.h"  // for Point

#include <vector>                // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class RealPin; }
namespace dealii { template <int dim, int spacedim> class Triangulation; }
#endif

namespace Forest
{
  using namespace dealii;

  /// @todo this name is provisional. Something better?
  template <int dim>
  class RealAssembly
  {
  public:

    /**
     * @brief Empty constructor
     */
    RealAssembly ();

    /**
     * @brief Used to initialize the object
     * @param lattice
     * @param origin
     * @param correction_factor
     */
    void
    init (const Forest::Lattice & lattice,
          Point<dim> origin,
          const double correction_factor);

    /**
     *
     * @param pins
     * @param real_pins
     */
    void
    build_data (const std::map<unsigned int, Pin> & pins,
                std::vector<RealPin<dim> > & real_pins);

    /**
     *
     * @param real_pins
     * @param tria
     */
    void
    build_mesh (std::vector<RealPin<dim> > & real_pins,
                Triangulation<dim, dim> & tria);

    /**
     * Printing the information of the particular assembly
     */
    void
    info () const;

    /**
     * @brief Set the correction factor to @p correction_factor
     * @todo Study if this is necessary, and remove it if not
     * @param correction_factor
     */
    inline void
    set_correction_factor (double correction_factor)
    {
      m_correction_factor = correction_factor;
    }

    /**
     * @brief Return the number of pins for this assembly
     * @return number of pins
     */
    inline unsigned int
    get_n_pins () const
    {
      return m_n_pins;
    }

    unsigned int m_empty_pin = static_cast<unsigned int> (-1);

  private:

    /**
     * @brief Used to fill the object with data from a Forest::Lattice
     * object
     * @param lattice
     */
    void
    set_lattice_data (const Forest::Lattice & lattice);

    /**
     * Set the origin of the assembly to @p origin
     * @param origin
     */
    inline void
    set_origin (const Point<dim> origin)
    {
      m_origin = origin;
    }

    /**
     *
     * @param pins
     * @param real_pins
     */
    void
    build_pins (const std::map<unsigned int, Pin> & pins,
                std::vector<RealPin<dim> > & real_pins);

    /**
     *
     * @param pins
     * @param real_pins
     */
    void
    set_assembly_pin_map (const std::map<unsigned int, Pin> & pins,
                          std::vector<RealPin<dim> > & real_pins);

    /**
     *
     * @param real_pins
     */
    void
    set_assembly_grid (std::vector<RealPin<dim> > & real_pins);

    /**
     *
     */
    void
    build_additional_data ();

    /// @todo I should change this to private
  public:

    std::vector<unsigned int> m_n_nodes;
    std::vector<std::vector<unsigned int> > m_components;
    Forest::LattType m_type;
    //Point<dim> m_pitch;
    std::vector<std::vector<double> > m_length;

    Point<dim> m_origin;

    unsigned int m_pins_ini;
    unsigned int m_pins_end;
    unsigned int m_n_pins;

    double m_correction_factor;

    bool m_initialized;
    bool m_ready_to_build_mesh;
  };
}
#endif
