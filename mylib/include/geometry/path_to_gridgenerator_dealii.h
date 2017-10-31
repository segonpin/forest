#ifndef path_to_gridgenerator_sebas_h
#define path_to_gridgenerator_sebas_h

#include "deal.II/base/point.h"

#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim> class Triangulation; }
#endif

namespace dealii
{
  namespace GridGenerator
  {

    void
    extrude_triangulation (const Triangulation<1, 1> &input,
                           const unsigned int n_slices,
                           const double height,
                           Triangulation<2, 2> &result);

    /**
     * @brief Function to flatten the triangulation
     * @note This function works as far as we do not have hanging nodes
     * @param in_tria
     * @param out_tria
     */
    template <int dim, int spacedim1 = dim, int spacedim2 = spacedim1>
    void
    flatten_triangulation (const Triangulation<dim, spacedim1> &in_tria,
                           Triangulation<dim, spacedim2> &out_tria);

    /**
     * @brief Overload of merge_triangulation to deal with a vector of
     * triangulations (not just two).
     * @param tria
     * @param result
     */
    template <int dim, int spacedim = dim>
    void
    merge_triangulations (
      const std::vector<Triangulation<dim, spacedim> > &tria,
      Triangulation<dim, spacedim> &result);

    /**
     * @brief Function to build the mesh for a pin of type PinType::PIN
     * @details We assume the pin cell is a square of
     * length @p pitch with a circular subdomain of radius
     * @p pin_radius, centered at @p center.
     * @param tria
     * @param center
     * @param pitch
     * @param pin_radius
     * @param mat_id
     */
    template <int dim>
    void
    pin_cell (Triangulation<dim, dim> & tria,
              const Point<dim> & center,
              const Point<dim> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id);

    /**
     * @brief Function to build the mesh for a pin of type PinType::PIN
     * @details We assume the pin cell is a square of
     * length @p pitch with a circular subdomain of radius
     * @p pin_radius, centered at @p center.
     * @param tria
     * @param center
     * @param pitch
     * @param pin_radius
     * @param mat_id
     */
    template <int dim>
    void
    pin_cell2 (Triangulation<dim, dim> & tria,
               const Point<dim> & center,
               const Point<dim> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id);

    /**
     * @brief Function to build the mesh for a pin of type PinType::BOX
     *
     * @details We assume the pin cell is a square of length @p pitch
     * centered at @p center.
     *
     * @param tria
     * @param center
     * @param pitch
     * @param mat_id
     */
    template <int dim>
    void
    pin_box (Triangulation<dim, dim> & tria,
             const Point<dim> & center,
             const Point<dim> & pitch,
             const std::vector<unsigned int> & mat_id);
  }
}

#endif
