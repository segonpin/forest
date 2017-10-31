/**
 * @brief  PostProcessing class template declarations
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   output/forest_post_processing.h
 */

#ifndef FOREST_POST_PROCESSING_H
#define FOREST_POST_PROCESSING_H

#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&
#include "neutronics/state_fwd.h"       // For class State&
#include "input/input_fwd.h"            // For class Input&
#include "deal.II/lac/vector.h"         // for Vector
#include "algebra/equation_fwd.h"       // For class Equation&

#include <string>                       // for string
#include <vector>                       // for vector

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
namespace dealii { template <int dim, int spacedim> class FiniteElement; }
#endif

/**
 * @brief Classes related to post-processing and printing the solution
 */
namespace Forest
{
  using namespace dealii;

  /**
   * @brief PostProcessing
   */
  template <int dim>
  class PostProcessing
  {
  public:

    /** @brief Constructor */
    PostProcessing (const State<dim> &state);

    /** @brief calculate the scaling factor to have average one. */
    double
    calculate_scaling(const std::vector<double> & vec);

    /** @brief Build the pin power distribution */
    void
    build_flux_power_data ();

    /** @brief Print the pin power distribution */
    void
    print_pin_data (std::string & bin,
                    const std::vector<double> & quantity,
                    const std::string &filename,
                    const unsigned int n_indent);

    /** @brief Print the pin power distribution */
    void
    print_ass_data (std::string & bin,
                    const std::vector<double> & quantity,
                    const std::string &filename,
                    const unsigned int n_indent);

    /** @brief Print the pin power distribution */
    void
    print_xml_data (const std::string &filename);

    /**
     * @brief Printing to the vtk file
     * @details We print the following
     *   - The fission rate associated with each dofs
     *   - The scalar flux associated with each dofs
     *   - The materials associated with each cell
     *   - The `user_index` associated with each cell
     */
    void
    print_vtk () const;

    /** @brief Printing to the eps file */
    void
    print_eps () const;

    /** @brief Printing to the gpl file */
    void
    print_gpl () const;

    /**
     * @brief Printing the output
     * @details Write the solution in eps for 1d, and to vtk for higher
     * dimensions (if the output flag is set to yes).
     * I also write the mesh in vtk for 2d and 3d.
     * @todo No mesh in 1d. In eps should be adding some extra dimension
     * (to plot like a bar) and in vtk in the same way?
     */
    void
    print_image_data () const;

    /**@brief Return the scaling factor to be used by other class */
    double
    get_scaling() const
    {
      return m_scaling;
    }

  private:

    const State<dim> & mp_state;
    const Input & data;
    const ProbGeom<dim> & geom;
    const FiniteElement<dim, dim> & fe;
    const DoFHandler<dim, dim> & dof_handler;


    std::vector<double> m_power_cell;
    std::vector<std::vector<double> > m_flux_cell;
    std::vector<double> m_power_pin;
    std::vector<std::vector<double> > m_flux_pin;
    std::vector<double> m_power_ass;
    std::vector<std::vector<double> > m_flux_ass;

    Vector<double> m_power_dofs;

    std::vector<double> m_weight_cell;
    std::vector<double> m_weight_pin;
    std::vector<double> m_weight_ass;

    double m_scaling;

    std::vector<double> m_norm_power_cell;
    std::vector<std::vector<double> > m_norm_flux_cell;
    std::vector<double> m_norm_power_pin;
    std::vector<std::vector<double> > m_norm_flux_pin;
    std::vector<double> m_norm_power_ass;
    std::vector<std::vector<double> > m_norm_flux_ass;

    Vector<double> m_norm_power_dofs;
    std::vector<Vector<double> > m_norm_flux_dofs;

    bool m_verbose;
  };

} /* end of namespace Forest */

#endif /* FOREST_POST_PROCESSING_H */
