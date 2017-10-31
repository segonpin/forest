/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input_geom.h
 * @brief  InputGeom class template declarations
 */
#ifndef FOREST_INPUT_GEOM_H
#define FOREST_INPUT_GEOM_H

#include <string>
#include <vector>
#include <map>

// This is for the xml parser from boost

namespace Forest
{
  /**
   * @brief Type of pin
   * @ingroup ForestInput
   */
  enum class PinType
  {
    BOX, /**< Rectangular box */
    PIN, /**<  Pin cell with circular pin fuel */
    PINNEW, /**< Pin cell with circular pin fuel */
    PINBOX /**< Pin cell with rectangular pin fuel */
  };

  /**
   * @brief Type of lattice
   * @ingroup ForestInput
   */
  enum class LattType
  {
    GRID, /**<  Lattice composed by Pin_type::box */
    PINMAP /**<  Lattice composed by Pin_type::pin and Pin_type::box */
  };

  /**
   * @class Core
   * @brief Container for the data of the core
   * @ingroup ForestInput
   */
  class Core
  {
  public:
    Core ()
    :m_max_dim(3)
    {};

    Core (const std::vector<unsigned int> & n_nodes,
         const std::vector<std::vector<double> > & length,
         const std::vector<std::vector<unsigned int>> & composition,
         const std::vector<unsigned int> & bcs,
         const std::string & name = "Core")
    {
      set(n_nodes, length, composition, bcs, name);
    };

    void
    set (const std::vector<unsigned int> & n_nodes,
         const std::vector<std::vector<double> > & length,
         const std::vector<std::vector<unsigned int>> & composition,
         const std::vector<unsigned int> & bcs,
         const std::string & name = "Core");

    /** @brief Maximum spatial dimension allowed. */
    unsigned int m_max_dim = 3;
    /** @brief Empty Assembly. */
    unsigned int m_empty_assembly = static_cast<unsigned int> (-1);
    /** @brief number of global refinements to be applied before starting */
    //unsigned int m_n_ref;

    /** @brief name for the core (not really used) */
    std::string m_name;

    /**
     * @brief Number of nodes per axis.
     * @details This variable has @p dim components, being used as:
     *  -# n_nodes[0] = number of nodes in x-axis
     *  -# n_nodes[1] = number of nodes in y-axis
     *  -# n_nodes[2] = number of nodes in z-axis
     * @note Not all the nodes in the lattice are necessarily not empty
     */
    std::vector<unsigned int> m_n_nodes;

    /**
     * @brief Number of nodes.
     * @details This variable has @p dim vectors, having vector @p d
     * size length[d].size() = n_nodes[d], being used as:
     *  -# length[0][i] = length of node @p i in x-axis
     *  -# length[1][i] = length of node @p i in y-axis
     *  -# length[2][i] = length of node @p i in z-axis
     * @note Not all the nodes in the lattice are necessarily not empty
     */
    std::vector<std::vector<double> > m_length;

    /**
     * @brief Each non-empty position has a unique index.
     * @details This variable has @p dim indices, being the first
     * index for the z-component, second for the y-component, and the third
     * one for the x-axis component:
     *  -# indices[k][:][:] = x-y plane at height @p k
     *  -# indices[k][j][:] = line parallel to x-axis at height @p k
     *  and ordinate @p j
     * @note z and y-components run downwards from n_nodes[2] and n_nodes[1]
     * until zero, respectively (only x-component increases).
     */
    std::vector<std::vector<std::vector<unsigned int> > > m_indices;

    /**
     * @brief the component for each assembly defined by its index
     */
    std::vector<unsigned int> m_components;

    /**
     * @brief boundary conditions
     */
    std::vector<unsigned int> m_bcs;

  private:

    /**
     * @brief assigning n_nodes (for dimensions higher than dim it is filled with 1's)
     * @param bin
     */
    void
    set_n_nodes (std::string & bin);

    void
    set_n_nodes (const std::vector<unsigned int> & n_nodes);

    void
    set_length (std::string & bin);

    void
    set_length (const std::vector<std::vector<double> > & length);

    void
    set_components (std::string & bin);

    void
    set_components (const std::vector<std::vector<unsigned int> > & composition);

    void
    set_bcs (std::string & bin);

    void
    set_bcs (const std::vector<unsigned int> & bcs);

  };

  /**
   * @class Lattice
   * @brief Container for the data of an assembly
   * @ingroup ForestInput
   */
  class Lattice
  {
  public:
    Lattice()
    : m_id(0),
      m_latt_type(LattType::GRID)
    {};

    Lattice(const unsigned int id,
            const std::vector<unsigned int> & n_nodes,
            const std::vector<std::vector<unsigned int> > & components,
            const std::vector<std::vector<double> > & length,
            const std::vector<double> & water_gap,
            const std::string & latt_type,
            const std::string & name = "Assembly")
    {
      set (id, n_nodes, components, length, water_gap, latt_type, name);
    };


    void
    set (const unsigned int id,
         const std::vector<unsigned int> & n_nodes,
         const std::vector<std::vector<unsigned int> > & components,
         const std::vector<std::vector<double> > & length,
         const std::vector<double> & water_gap,
         const std::string & latt_type,
         const std::string & name = "Assembly");


    /** @brief Maximum spatial dimension allowed. */
    static const unsigned int m_max_dim = 3;
    /** @brief Empty pin. */
    unsigned int m_empty_pin = static_cast<unsigned int> (-1);
    /** @brief identification for this lattice */
    unsigned int m_id;
    /** @brief type of lattice. */
    LattType m_latt_type;
    /** @brief name of the lattice (purely informative). */
    std::string m_name;

    /**
     * @brief Number of nodes per axis.
     * @details This variable has @p dim components, being used as:
     *  -# n_nodes[0] = number of nodes in x-axis
     *  -# n_nodes[1] = number of nodes in y-axis
     *  -# n_nodes[2] = number of nodes in z-axis
     * @note Not all the nodes are necessarily not empty
     */
    std::vector<unsigned int> m_n_nodes;

    /** @brief the component for each pin defined by its index */
    std::vector<std::vector<unsigned int> > m_components;

    /** @brief Minimum refinement for this assembly */
    //unsigned int m_n_ref;

    /**
     * @brief length of nodes (for Latt_type::grid)
     * @details This variable has @p dim vectors, having vector @p d
     * size length[d].size() = n_nodes[d], being used as:
     *  -# length[0][i] = length of node @p i in x-axis
     *  -# length[1][i] = length of node @p i in y-axis
     *  -# length[2][i] = length of node @p i in z-axis
     * @note Not all the nodes in the lattice are necessarily not empty
     */
    std::vector<std::vector<double> > m_length;

    /** @brief pitch and water_gap are used with  Latt_type::pin_map */
    //std::vector<double> m_pitch;
    /** @brief water gap.*/
    std::vector<double> m_water_gap;

    /** @brief vector with the size (sum of all pin sizes) of the assembly */
    std::vector<double> m_size;

  private:

    /**
     * @brief assigning n_nodes
     * @details It is always a full size vector (max_dim size),
     * so for dimensions higher than dim it is filled with 1's .
     */
    void
    set_n_nodes (const std::string & bin);

    void
    set_n_nodes (const std::vector<unsigned int> & n_nodes);

    void
    set_length (const std::string & bin);

    void
    set_length (const std::vector<std::vector<double> > & length);

    void
    set_components (const std::string & bin);

    void
    set_components (const std::vector<std::vector<unsigned int> > & components);

    void
    set_type (const std::string & bin);

    LattType
    parse (const std::string & bin);

    void
    set_size_from_length ();

  };

  /**
   * @class Pin
   * @brief Container for the data of a single pin
   * @ingroup ForestInput
   */
  class Pin
  {
  public:
    Pin()
    :m_id(0),
     m_pin_type(PinType::BOX),
     m_fuel_radius(0.0)
    {};

    Pin(const unsigned int id,
         const std::vector<unsigned int> & materials,
         const std::string & pin_type,
         const double fuel_radius,
         const std::string & name = "Pin")
    {
      set(id, materials, pin_type, fuel_radius, name);
    }

    void
    set (const unsigned int id,
         const std::vector<unsigned int> & materials,
         const std::string & pin_type,
         const double fuel_radius,
         const std::string & name = "Pin");

  //private:
    /** @brief Minimum refinement for this pin */
    //unsigned int m_min_ref;

    /** @brief name for the pin (not really used) */
    std::string m_name;

    /** @brief identification for this lattice */
    unsigned int m_id;

    /** @brief type of pin. */
    PinType m_pin_type;

    /** @brief radius for the pin when a circular pin is considered */
    double m_fuel_radius;

    /** @brief materials composing the pin
     * @details materials is a vector of two components, being the first
     * one the material_id for the water, and the second one the
     * material_id for the fuel
     */
    std::vector<unsigned int> m_materials;

  private:

    PinType
    parse (const std::string & bin);

    void
    set_type (std::string & bin);

  };

  /**
   * @class InputGeom
   * @brief Class containing input data to construct the geometry
   * @ingroup ForestInput
   */
  class InputGeom
  {
  public:
    /** Constructor */
    InputGeom ();
    InputGeom (const std::string &filename);
    /** @brief spatial dimension for the problem */
    unsigned int m_dim;
    /** @brief Container for the information of the core */
    Core m_core;
    /** @brief Container for the information of the different lattices */
    std::map<unsigned int, Lattice> m_lattices;
    /** @brief Container for the information of the different pins */
    std::map<unsigned int, Pin> m_pins;

    void
    set (unsigned int dim,
         Core core,
         std::map<unsigned int, Lattice> lattices,
         std::map<unsigned int, Pin> pins);

    unsigned int
    get_dim () const
    {
      return m_dim;
    }

    /**
     * @brief read and store data
     * @details We read the file @p filename from disk and store the
     * information in the class, using the different
     * structures provided.
     * @param[in] filename
     */
    void
    load_xml (const std::string &filename);

    /*void
    set_lattice_grid (Lattice & lattice);*/

    InputGeom
    generate_from_lattice(const unsigned int ass_id) ;

    /**
     * @brief write data
     * @details We write the file @p filename in the disk from the
     * information provided in the different structures provided.
     * @param[in] filename
     */
    void
    save_xml (const std::string &filename);

    /** @brief Check the different structures */
    std::string
    show ();

    /** @brief Check the different structures */
    void
    check () const;

    /** @brief Check the core */
    void
    check_core () const;

    /** @brief Check the lattices */
    void
    check_lattices () const;

  };

} // end of namespace Forest

#endif /* FOREST_INPUT_GEOM_H */
