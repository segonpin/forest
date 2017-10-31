/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/extract_bv.h
 * @brief  ExtractBV class template declarations
 */

#ifndef FOREST_EXTRACT_BV_H
#define FOREST_EXTRACT_BV_H

#include "neutronics/state_fwd.h"
#include "algebra/sn_vector.h"          // for SnVector
#include "angle/quadrature_base_fwd.h"
#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&

#include "deal.II/dofs/dof_handler.h"   // for DoFHandler

#include <map>                          // for map
#include <set>                          // for set
#include <memory>                       // for shared_ptr
#include <string>                       // for string
#include <vector>                       // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class Homogenization; }
namespace Forest { template <int dim> class HomData; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class ExtractBV
   * @ingroup ForestNeutronics
   * Extracts boundary values for the system.
   */
  template <int dim>
  class ExtractBV
  {
  public:
    /**
     @brief Constructor.
     @details This is not an empty constructor, but only initialize references.
     It should be initialize before the object is used for something else.
     @param state
     */
    ExtractBV (const State<dim> &state);

    /**
     @brief Intialize the members of the class that need it.
     */
    void
    initialize ();

  private:

    /** @brief Indicating the level of verbose: 0 = none; 1 = everything. */
    bool m_verbose;

    /** @brief Reference to the state object */
    const State<dim> &mp_state;

    /** @brief Angular quadrature */
    std::shared_ptr<QuadratureBase<dim> > m_ord;

    /** @brief Geometry */
    const ProbGeom<dim> &mp_geom;
    const DoFHandler<dim, dim> &m_dof_handler;
    const unsigned int m_n_groups;
    const unsigned int m_n_angles;

    /// @todo move this value to the materials input?
    static const unsigned int m_n_leg_mom = 1;

  public:

    /**
     @brief Return the memory consumption of an object of this class.
     @return The memory consumption
     @todo prepare this function for the members of this class?
     */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      return memory_consumption;
    }

    /**
     * @brief Identifying and storing the boundary cells for each pin.
     * @details The information is store in m_pin_bc,
     * which is a map containing the pairs
     * @code{.cpp}
     * std::pair<c_iter, unsigned int> (cell,f)
     * @endcode
     * for easy access through an iterator
     * @code{.cpp}
     * typename std::map<c_iter, unsigned int>::iterator it, end;
     * it = m_pin_bc[pin][face].begin ();
     * end = m_pin_bc[pin][face].end ();
     * for (; it != end; ++it)
     * {
     *   c_iter cell = it.first;
     *   unsigned int face = it.second;
     * }
     * @endcode
     */
    void
    set_pin_bc ();

    /** @brief Fill the boundary weights. */
    void
    set_pin_bw ();

    /**
     * @brief Extract the boundary conditions for the pins.
     * @details We take the scalar flux and its moments, we converted to the
     * angular flux with the moments-to-directions operator, and extract the
     * flux from there.
     * @warning This should be done when the flux is already converged, so the
     * angular flux has the right value at the boundaries that can be extracted,
     * without solving again the problem. This is based on the assumption that
     * \f{align}{
     * MD \approx I
     * \f}
     * so we have the approximation
     * \f{align}{
     * \Psi\approx MD \Psi = M \phi ,
     * \f}
     * what in general is absolutely ** not ** true.
     */
    void
    set_pin_bv ();

    /** @brief Fill the containers for the bcs. with its values. */
    void
    set_pin_bv_g_n (const Vector<double> & psi_n,
                    const unsigned int g,
                    const unsigned int n);

    /**
     * @brief Identifying and storing the boundary cells for each assembly.
     * @details The information is store in m_ass_bc,
     * which is a map containing the pairs
     * @code{.cpp}
     * std::pair<c_iter, unsigned int> (cell,f)
     * @endcode
     * for easy access through an iterator
     * @code{.cpp}
     * typename std::map<c_iter, unsigned int>::iterator it, end;
     * it = m_ass_bc[ass][face].begin ();
     * end = m_ass_bc[ass][face].end ();
     * for (; it != end; ++it)
     * {
     *   c_iter cell = it.first;
     *   unsigned int face = it.second;
     * }
     * @endcode
     */
    void
    set_ass_bc ();

    /**
     * @brief Identifying the pin indices over the boundary of the assembly.
     * @details The information is store in m_ass2pin,
     * which is a map containing the pairs
     * @code{.cpp}
     * std::pair<unsigned int, unsigned int> (pin, counter)
     * @endcode
     * for easy access through an iterator
     * @code{.cpp}
     * typename std::map<unsigned int, unsigned int>::iterator it, end;
     * it = m_ass2pin[ass][face].begin ();
     * end = m_ass2pin[ass][face].end ();
     * for (; it != end; ++it)
     * {
     *   unsigned int pin = it.first;
     *   unsigned int counter = it.second;
     * }
     * @endcode
     */
    void
    set_ass2pin ();

    /** @brief print assembly boundary pointers. */
    void
    print_ass_bc ();

    /** @brief Fill the boundary weights. */
    void
    set_ass_bw ();

    /** @brief Fill the boundary values. */
    void
    set_ass_bv ();



    /**
     * @brief Calculate neutron current
     * @details Knowing the angular flux \f$ \Psi \f$, the net current is
     * readily available through integration as follows
     *
     * \f{align*}{
     * \vec{J}
     * & = \int \Omega \Psi(\mathbf{r},\Omega) \ \text{d}\Omega
     *  \approx \sum_{n} \omega_{n} \Omega_{n} \Psi_{n}(\mathbf{r})
     * \f}
     *
     * In a similar way, the partial currents for incoming and outgoing
     * directions over a face can be calculated as
     *
     * \f{align*}{
     * J_{-}
     * & = \int_{\vec{n} \cdot \Omega  < 0}
     * (- \vec{n} \cdot \Omega)\Psi(\mathbf{r},\Omega) \ \text{d}\Omega
     * \approx \sum_{\vec{n} \cdot \Omega  < 0}
     * \omega_{n} (- \vec{n} \cdot \Omega_{n}) \Psi_{n}(\mathbf{r})
     * \f}
     *
     * and
     *
     * \f{align*}{
     * J_{+}
     * & = \int_{\vec{n} \cdot \Omega  > 0}
     * (\vec{n} \cdot \Omega) \Psi(\mathbf{r},\Omega) \ \text{d}\Omega
     * \approx \sum_{\vec{n} \cdot \Omega  > 0}
     * \omega_{n} (\vec{n} \cdot \Omega_{n}) \Psi_{n}(\mathbf{r})
     * \f}
     *
     * where \f$ \vec{n} \f$ is the outward vector over the boundary of the
     * domain. We ca use the previous definitions to calculate the normal
     * component of the net current over the direction of the outward direction
     * perpendicular to the face as follows
     *
     * \f{align*}{
     * \vec{J} \cdot \vec{n}
     * & = \int (\vec{n}\cdot \Omega)\Psi(\mathbf{r},\Omega) \ \text{d}\Omega \\
     * & = \int_{\vec{n} \cdot \Omega  > 0}
     * (\vec{n} \cdot \Omega) \Psi(\mathbf{r},\Omega) \ \text{d}\Omega
     * - \int_{\vec{n} \cdot \Omega  < 0}
     * (- \vec{n} \cdot \Omega) \Psi(\mathbf{r},\Omega) \ \text{d}\Omega
     * = J^{+} - J^{-}
     * \f}
     *
     * @todo The angular flux is obtained approximately from the scalar flux.
     * It should be correctly implemented with the right angular flux.
     */
    void
    calculate_current_pin ();

    void
    calculate_current_ass ();

    void
    calculate_moments_pin ();

    void
    calculate_moments_ass ();

    /** @brief Generate the coarser angular representation for the bv */
    void
    coarsening_angular ();

    /**
     * @brief Generate the coarser spatial representation for the bv
     * @todo Automate the coarsening level and make it generic
     */
    void
    coarsening_spatial ();

  public:

    /**
     * @brief read the bv from a file.
     * @todo Change BOOST_FOREACH for new ranged-based for?
     * ptree d = pt.get_child("boundary_conditions");
     * //for (const auto &face : d)
     * for (const std::pair<std::string, ptree> &face : d)
     * {
     * if (face.first != "face")
     * continue;
     * }
     */
    void
    read_ass_bv ();

    /**
     * @brief print the bv to a file.
     * @todo Change BOOST_FOREACH for new ranged-based for?
     * ptree d = pt.get_child("boundary_conditions");
     * //for (const auto &face : d)
     * for (const std::pair<std::string, ptree> &face : d)
     * {
     * if (face.first != "face")
     * continue;
     * }
     */
    void
    print_ass_bv (const std::string &filename_prefix) const;

  protected:

    /** @brief shortcut for std::vector<double> */
    typedef std::vector<double> stdvd1;
    /** @brief recursive shortcut for std::vector<std_vd1> */
    typedef std::vector<stdvd1> stdvd2;
    /** @brief recursive shortcut for std::vector<std_vd2> */
    typedef std::vector<stdvd2> stdvd3;
    /** @brief recursive shortcut for std::vector<std_vd3> */
    typedef std::vector<stdvd3> stdvd4;

    /** @brief shortcut for std::vector<unsigned int> */
    typedef std::vector<unsigned int> std_vui1;
    /** @brief recursive shortcut for std::vector<std_vui1> */
    typedef std::vector<std_vui1> std_vui2;
    /** @brief recursive shortcut for std::vector<std_vui2> */
    typedef std::vector<std_vui2> std_vui3;
    /** @brief recursive shortcut for std::vector<std_vui3> */
    typedef std::vector<std_vui3> std_vui4;

    /** @brief shortcut for the cell_iterator type */
    typedef typename DoFHandler<dim, dim>::active_cell_iterator c_iter;

    /** @brief shortcut for the map <cell_iterator, face>  */
    using map_cf = std::map<c_iter, unsigned int>;

    /**
     * m_pin_bc[pin][face] is a map containing the pairs
     * pair(cell,f), so we can easily find the cells and faces inside a
     * triangulation defining the boundary @p face of the pin @p pin
     */
    std::vector<std::vector<map_cf> > m_pin_bc;

    /**
     * @brief Container for the pin faces boundary weight
     * @details The arguments here are the @p pin_index and
     * the @p pin_face, storing the measure of the face as follows
     * @code
     * double boundary_weight = m_pin_bw[pin][f];
     * @endcode
     */
    stdvd2 m_pin_bw;

    /**
     * @brief Container for the pin faces boundary values
     * @details The arguments here are the @p pin_index,
     * the @p pin_face, the @p energy_group and the @p angle.
     * For each of this combinations, we have the average value of the flux at
     * the specified face of the pin as follows
     * @code
     * double boundary_value = m_pin_bv[pin][f][group_g][angle];
     * @endcode
     */
    stdvd4 m_pin_bv;

    /**
     * @brief Container for the pin faces boundary values
     * @details The arguments here are the @p pin_index,
     * the @p pin_face, the @p energy_group and the @p moment.
     * For each of this combinations, we have the average value of the flux at
     * the specified face of the pin as follows
     * @code
     * double boundary_value_m = m_pin_bvm[pin][f][group_g][moment];
     * @endcode
     */
    stdvd4 m_pin_bvm;

    /**
     * @brief Pins over the boundary of the assembly
     * @details assebmly_boundary[ass_i][face_i] is a set containing the
     * cells, so we can easily find the cells and faces inside a
     * triangulation defining the boundary @p face_i of the assembly @p ass_i
     * @note used in inner_iteration and boundaryvalues
     */
    std::vector<std::vector<std::set<unsigned int> > > m_ass2pin;

    /**
     * @brief boundary container for the assembly sorting cells and faces
     * @details m_ass_bc[ass_i][face_i] is a map containing the pairs
     * pair(cell,f), so we can easily find the cells and faces inside a
     * triangulation defining the boundary @p face_i of the assembly @p ass_i
     */
    std::vector<std::vector<map_cf> > m_ass_bc;

    /**
     * @brief Container for the ass faces boundary weight.
     * @details
     * @code
     * double boundary_weight = m_ass_bw[ass][f];
     * @endcode
     */
    stdvd2 m_ass_bw;

    /**
     * @brief Container for the boundary values for the assembly.
     * @details
     * @code
     * double bv = m_ass_bv[ass][f][group_g][angle];
     * @endcode
     */
    stdvd4 m_ass_bv;

    /**
     * @brief Container for the boundary values for the assembly.
     * @details
     * @code
     * double bvm = m_ass_bvm[ass][f][group_g][moment];
     * @endcode
     */
    stdvd4 m_ass_bvm;

    /**
     * @brief Container for the boundary values for the flux.
     * @details
     * @code
     * double bv = m_phi0[pin][f][group_g];
     * @endcode
     */
    stdvd3 m_phi0_pin;
    stdvd3 m_phi0_ass;

    /** @brief Container for the values for the outgoing partial current.
     * @details
     * @code
     * double bv = m_jp_o[pin][f][group_g];
     * @endcode
     */
    stdvd3 m_jp_o_pin;
    stdvd3 m_jp_o_ass;

    /** @brief Container for the values for the incoming partial current.
     * @details
     * @code
     * double bv = m_jp_i[pin][f][group_g];
     * @endcode
     */
    stdvd3 m_jp_i_pin;
    stdvd3 m_jp_i_ass;

    /** @brief Container for the boundary values for the net current.
     * @details
     * @code
     * double bv = m_jnet[pin][f][group_g][dim];
     * @endcode
     */
    stdvd4 m_jnet_pin;
    stdvd4 m_jnet_ass;

    /** @brief Container for the boundary values for the net current.
     * @details
     * @code
     * double bv = m_j[pin][f][group_g];
     * @endcode
     */
    stdvd3 m_j_pin;
    stdvd3 m_j_ass;


    /** @brief Container for the boundary values for the second moment flux.
     * @details
     * @code
     * double bv = m_phi2[pin][f][group_g];
     * @endcode
     */
    stdvd3 m_phi2_pin;
    stdvd3 m_phi2_ass;

  public:

    /** @brief Contains the spherical harmonics order for the bcs. */
    std::vector<unsigned int> m_bcs_order;

    /** @brief Scale the interface values with a given scaling factor*/
    void
    scaling(const double scaling);
    /** @brief Scale the interface values with a given scaling factor*/
    void
    scaling_moments(const double scaling);

    friend class Homogenization<dim>;
    friend class HomData<dim>;

  };
} // end of namespace Forest

#endif /* FOREST_EXTRACT_BV_H_ */
