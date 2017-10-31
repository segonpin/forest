#ifndef REAL_PIN_H
#define REAL_PIN_H

#include "input/input_geom.h"

#include "deal.II/base/point.h"

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim> class Triangulation; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @brief Contains the information to build the triangulation for a pin
   */
  template <int dim>
  class RealPin
  {
  public:
    /**
     * @brief Empty constructor
     */
    RealPin ();

    /**
     * @brief The initializer
     * @details This initialized is provided we have a pin object
     * @param pin
     * @param correction_factor
     */
    void
    init (const Forest::Pin & pin,
          const double correction_factor = 1.0);

    /**
     * @brief The initializer
     * @details Here we do not have a pin object (like in a lattice_type::grid)
     * @param material_id
     */
    void
    init (const unsigned int material_id);

    /**
     * @brief set the pitch
     * @param pitch
     */
    inline void
    set_pitch (const Point<dim> pitch)
    {
      m_pitch = pitch;
    }

    /**
     * @brief Set the center of the pin
     * @param pin_origin
     */
    inline void
    set_center (const Point<dim> pin_origin)
    {
      m_center = pin_origin + 0.5 * m_pitch;
    }

    /**
     * @brief Get the center of the pin
     * @return
     */
    inline Point<dim>
    get_center () const
    {
      return m_center;
    }

    /**
     * @brief Get the center of the pin
     * @return
     */
    inline Point<dim>
    get_pitch () const
    {
      return m_pitch;
    }

    /**
     * @brief Get the radius
     */
    inline double
    get_radius () const
    {
      return m_radius;
    }

    /**
     * @brief Return the type of the pin
     * @return
     */
    inline Forest::PinType
    get_type () const
    {
      return m_type;
    }

    /**
     * @brief Check whether the point @p p is inside the pin or not
     * @param p
     * @return
     */
    bool
    is_inside (const Point<dim> p) const;

    /**
     * @brief Check the infinite norm from @p point to the center of the pin in xy
     * @param point
     * @return
     */
    double
    norm_inf_xy (const Point<dim> & point) const;

    /**
     * @brief Build a Triangulation<dim> object for this pin
     * @details Depending on the type of pin, different generators are used. It
     * should be possible to easily extend this definition to more type of pins.
     * @note The triangulations associated with the pins will then be merged
     * into a bigger triangulation (meaning that it represents an assembly or a
     * core). Thus, one basic requirement of this way of building the mesh is
     * that the pins can be merged together without having hanging nodes
     * (i.e. not possible to use non-matching meshes).
     * @param tria
     */
    void
    build_mesh (Triangulation<dim, dim> & tria);

    /**
     * Displays the information about this pin
     */
    void
    info () const;

  private:

    /// @p m_pitch  is the size of the pin when considered as a box.
    Point<dim> m_pitch;
    /// @p m_center denotes the center of the pin
    Point<dim> m_center;
    /// @p m_type is used to define the way the mesh will be constructed
    Forest::PinType m_type;
    /// @p m_materias is a vector of materials, the first one assigned to the
    /// fuel and the second one to the moderator.
    std::vector<unsigned int> m_materials;
    /// when the pin cell is of type @p Pin_type::pin, @p m_radius is the
    /// radius of the pin.
    double m_radius;
  };
}

#endif
