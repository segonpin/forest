#ifndef PROB_GEOM_H
#define PROB_GEOM_H

#include "input/input_fwd.h"            // For class Input&

#include "deal.II/grid/tria.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <memory>
#include <string>
#include <vector>

#ifndef DOXYGEN
namespace Forest { class Core; template <int dim> class RealPin; }
#endif

namespace Forest
{
  using namespace dealii;

  template <int dim>
  class ProbGeom
  {
  public:
    /**
     */
    ProbGeom (const Input & data);

    const Input & mp_data;

    std::vector<std::shared_ptr<Manifold<dim, dim> > > m_manifold;
    Triangulation<dim, dim> m_tria;

    // we use this to set the power map after the calculations
    std::vector<std::vector<unsigned int> > map_pin_index;

    /**
     * @brief this is defined identify an empty place (no pin or assembly).
     *
     * @details for example, the pin distribution can be defined by
     * @code
     * const std::vector<std::vector<unsigned int> > pins_map =
     *  {{ 1,     1,     1,     1    },
     *   { 1,     1,     1,     empty},
     *   { 1,     1,     empty, empty},
     *   { 1,     empty, empty, empty}};
     * @endcode
     */
    const Core & m_core;


    /** @brief We print the mesh */
    void
    print_mesh_eps (const std::string & output) const;

    /** @brief We print the mesh */
    void
    print_mesh_vtk (const std::string & output) const;

  private:

    void
    set_boundary_indicators ();

    void
    set_map_pin_index ();

    /**
     * @brief It calculates a correction factor such that the area of the
     * pin with radius @p r approximated by a polygon defined by @p n_ref
     * is equal to the area of the circle.
     * @details The area of the circle is
     * \f[ A_c = \pi * R^2 \f]
     * and the area of the polygon after @p n_ref refinements is
     * \f[
     *   A_p = \frac{1}{2} * 4 *  2^{n_\text{ref}}
     *   * (R*\text{correction\_factor})^2 * sin (\frac{2*pi}{4*2^{n_\text{ref}}})
     * \f]
     * doing \f$ A_c = A_p \f$, and isolating @p correction_factor we obtain
     * \f[
     *    \text{correction\_factor}
     *    = \sqrt{\frac{\pi*2}{4* 2^{n_\text{ref}}}}*sin(\frac{2*pi}{4*2^{n_\text{ref}}})
     * \f]
     * simplifying we obtain
     * \f[
     *  \text{correction\_factor}
     *   = \sqrt{\frac{pi}{2^{n_\text{ref+1}}}*sin(pi/(2^{n_\text{ref+1}}))}
     * \f]
     * \f[
     *  \text{correction\_factor}
     *   = \sqrt(pi/(std::pow(2,n_\text{ref+1})*sin(pi/(std::pow(2,n_\text{ref+1})))))
     * \f]
     * @todo The formulas in the documentation are wrong. Must be fixed.
     */
    double
    area_correction_factor (const unsigned int n_ref);

    void
    set_n_assemblies ();

    void
    set_pins_per_assembly ();

    void
    set_n_pins ();

    /**
     * @brief Attach the manifold to the mesh and refine it.
     * @details We define the manifolds, associate them with the manifold_id
     * used in the mesh for the particular manifold (circular pin shape),
     * and refine the mesh.
     *
     * @bug For the manifold_id we will start with the number
     * @p first_manifold = 1024 (2^10). There is no reason to choose this number
     * in particular, but using a number smaller than 256 gives problems with
     * manifold_id == 256, so we skip the problem at the moment,
     * and we have to prepare an example to submit to dealii user groups.
     * This bug seems to be related with the use of unsigned char somewhere...
     */
    void
    set_manifold (Triangulation<dim, dim> & tria,
                  std::vector<RealPin<dim> > & real_pins,
                  const unsigned int first_manifold = 1024);

    /**
     * @todo rethink about how do I want to control the refinements to
     * allow local refinements when I have a pin with the circle, but allow
     * also refinements when the pin is a square box generated with the grid
     * option for the assembly (because now it is doing no refinement for this
     * case so I will add minimun global refinemens). */
    void
    refine_assembly (Triangulation<dim, dim> & tria,
                     std::vector<RealPin<dim> > & real_pins,
                     const unsigned int n_ref = 0);
    /**
     * @brief We set the different manifold_id's for the cells in different pins
     * @details given a triangulation @p tria, we run over all the cells of the
     * triangulation and compare their center with the center of the pins in order
     * to know if the particular cell represents the circular pin inside a pin
     * cell.
     * @note We do this after using @p extrude_triangulation because this
     * function does not preserve this information yet.
     * @note This is done in an inefficient way. We keep this algorithm for the
     * moment, but the memory consumption should be profiled to see if it needs
     * improvement.
     * @note Using cell-set_all_manifold_ids(manifold_ind) should work, but it
     * does not.
     * @param tria
     * @param real_pins
     * @param first_manifold
     */
    void
    assign_pin_manifolds (Triangulation<dim, dim> & tria,
                          std::vector<RealPin<dim> > & real_pins,
                          unsigned int first_manifold = 10);

    /**
     * @brief assign user_index (per pin) for the output of the power per pin
     * @details run over all the cells to assign the user_index (per pin)
     * for the output of the power per pin
     * @todo I have to add capabilities for the 3d
     * @todo it should be optimized (different search algorithm)
     */
    std::vector<unsigned int>
    assign_cell_index (Triangulation<dim, dim> & tria,
                       std::vector<RealPin<dim> > & real_pins);

    unsigned int n_assemblies;
    unsigned int n_pins;

    std::vector<unsigned int> pins_per_assembly;

    /**
     * @details pin_id = cells_to_pin[user_index] associate to the user_index
     * of a cell the pin_id for which they belong to.
     * @note the size for this vector should be the number of cells
     */
    std::vector<unsigned int> cell2pin;

    /**
     * @details assembly_id = pins_to_assembly[pin_id] associate to the
     * identifier of a pin the identifier of the assembly they belong to.
     * @note the size for this vector should be the number of pins
     */
    std::vector<unsigned int> pin2assembly;

    unsigned int m_mesh_regularity;

    void
    set_mesh_regularity ();


  public:

    unsigned int
    get_n_pins () const
    {
      return n_pins;
    }

    unsigned int
    get_n_assemblies () const
    {
      return n_assemblies;
    }

    unsigned int
    get_cell2pin (const unsigned int cell) const
    {
      if (cell == (unsigned int) -1)
        return -1;
      return cell2pin[cell];
    }

    unsigned int
    get_pin2assembly (const unsigned int pin) const
    {
      Assert(pin < n_pins, ExcMessage("Pin out of range"));
      return pin2assembly[pin];
    }

    unsigned int
    get_cell2assembly (const unsigned int cell) const
    {
      if (cell == (unsigned int) -1)
        return -1;
      return pin2assembly[cell2pin[cell]];
    }

    double
    memory_consumption () const;

    unsigned int
    get_mesh_regularity () const
    {
      return m_mesh_regularity;
    };
  };

} // end of namespace Forest

#endif
