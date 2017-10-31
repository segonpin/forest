/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input_geom.cc
 * @brief  Implementation of class template InputGeom
 */

#include <input/input_geom.h>

#include "utils/forest_utils_base.h"

#include <deal.II/base/exceptions.h>            // for Assert and ExcMessage

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION < 105500)
typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#else
typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#endif
//#include <boost/foreach.hpp>

#include <iostream>

namespace Forest
{

  void
  Core::set (const std::vector<unsigned int> & n_nodes,
             const std::vector<std::vector<double> > & length,
             const std::vector<std::vector<unsigned int>> & composition,
             const std::vector<unsigned int> & bcs,
             const std::string & name)
  {
    // core --------------------------------------------------------------
    set_n_nodes(n_nodes);
    set_length(length);
    set_components(composition);
    set_bcs(bcs);
    m_name = name;
  }

  void
  Core::set_n_nodes (std::string & bin)
  {
    std::vector<unsigned int> tmp_n_nodes;
    Forest::string_to_vector (bin, tmp_n_nodes);
    this->set_n_nodes(tmp_n_nodes);
  }

  void
  Core::set_n_nodes (const std::vector<unsigned int> & n_nodes)
  {
    this->m_n_nodes.resize (this->m_max_dim, 1);
    for (unsigned int i = 0; i < n_nodes.size (); ++i)
      this->m_n_nodes[i] = n_nodes[i];
  }

  /// @todo Here I should add some checks
  void
  Core::set_length (std::string & bin)
  {
    Forest::string_to_vector (bin, this->m_length);
  }

  /// @todo Here I should add some checks
  void
  Core::set_length (const std::vector<std::vector<double> > & length)
  {
    m_length = length;
  }

  /// @todo Here I should add some checks
  void
  Core::set_components (std::string & bin)
  {
    std::vector<std::vector<unsigned int> > composition;
    Forest::string_to_vector (bin, composition);
    set_components(composition);
  }

  void
  Core::set_components (const std::vector<std::vector<unsigned int> > & composition)
  {
    this->m_indices.push_back (composition);
    unsigned int index = 0;
    for (unsigned int j = 0; j < this->m_n_nodes[1]; ++j)
      for (unsigned int i = 0; i < this->m_n_nodes[0]; ++i)
        if (composition[j][i] != this->m_empty_assembly)
        {
          this->m_indices[0][j][i] = index;
          this->m_components.push_back (composition[j][i]);
          ++index;
        }
  }

  /// @todo Here I should add some checks
  void
  Core::set_bcs (std::string & bin)
  {
    Forest::string_to_vector (bin, this->m_bcs);
  }

  /// @todo Here I should add some checks
  void
  Core::set_bcs (const std::vector<unsigned int> & bcs)
  {
    m_bcs = bcs;
  }

  //-------------------------------------------------------------

  void
  Lattice::set (const unsigned int id,
                const std::vector<unsigned int> & n_nodes,
                const std::vector<std::vector<unsigned int> > & components,
                const std::vector<std::vector<double> > & length,
                const std::vector<double> & water_gap,
                const std::string & latt_type,
                const std::string & name)
  {
    m_id = id;
    set_n_nodes (n_nodes);
    set_components (components);
    set_length (length);
    set_type(latt_type);
    m_water_gap = water_gap;
    m_name = name;

    set_size_from_length ();
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_n_nodes (const std::string & bin)
  {
    std::vector<unsigned int> tmp_n_nodes;
    Forest::string_to_vector (bin, tmp_n_nodes);
    this->set_n_nodes(tmp_n_nodes);
  }

  void
  Lattice::set_n_nodes (const std::vector<unsigned int> & n_nodes)
  {
    this->m_n_nodes.resize (this->m_max_dim, 1);
    for (unsigned int i = 0; i < n_nodes.size (); ++i)
      this->m_n_nodes[i] = n_nodes[i];
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_length (const std::string & bin)
  {
    Forest::string_to_vector (bin, this->m_length);
  }

  void
  Lattice::set_length (const std::vector<std::vector<double> > & length)
  {
    m_length = length;
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_components (const std::string & bin)
  {
    Forest::string_to_vector (bin, this->m_components);
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_components (const std::vector<std::vector<unsigned int> > & components)
  {
    m_components = components;
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_type (const std::string & bin)
  {
     m_latt_type = this->parse(bin);
  }

  LattType
  Lattice::parse (const std::string & bin)
  {
     std::map<std::string, LattType> string2type;
     string2type["pin_map"] = LattType::PINMAP;
     string2type["grid"] = LattType::GRID;
     return string2type[bin];
  }

  /// @todo Here I should add some checks
  void
  Lattice::set_size_from_length ()
  {
    for (unsigned int i = 0; i < this->m_length.size (); ++i)
    {
      double total_length = 0;
      for (unsigned int j = 0; j < this->m_length[i].size (); ++j)
        total_length += this->m_length[i][j];
      this->m_size.push_back (total_length);
    }
  }

  /// @todo Here I should add some check
  void
  Pin::set_type (std::string & bin)
  {
    this->m_pin_type = this->parse(bin);
  }

  /// @todo Here I should add some check
  PinType Pin::parse (const std::string & bin)
  {
     std::map<std::string, PinType> string2type;
     string2type["pin"] = PinType::PIN;
     string2type["pinnew"] = PinType::PINNEW;
     string2type["box"] = PinType::BOX;
     string2type["pinbox"] = PinType::PINBOX;
     return string2type[bin];
  }

  void
  Pin::set (const unsigned int id,
            const std::vector<unsigned int> & materials,
            const std::string & pin_type,
            const double fuel_radius,
            const std::string & name)
  {
    m_name = name;
    m_id = id;
    m_pin_type = this->parse(pin_type);
    m_fuel_radius = fuel_radius;
    m_materials = materials;
  }

  InputGeom::InputGeom ()
      : m_dim (0), m_core(), m_lattices(), m_pins()
  {
  }

  InputGeom::InputGeom (const std::string &filename)
      : m_dim (0)
  {
    this->load_xml (filename);
    this->check ();
    // this->save_xml(settings.get_full_path_input()+"_out.geom.xml");
  }

  void
  InputGeom::set (unsigned int dim,
                  Core core,
                  std::map<unsigned int, Lattice> lattices,
                  std::map<unsigned int, Pin> pins)
  {
    m_dim = dim;
    m_core = core;
    m_lattices = lattices;
    m_pins = pins;
  }

  void
  InputGeom::load_xml (const std::string &filename)
  {
    unsigned int tmp_dim ;
    Core tmp_core;
    std::map<unsigned int, Lattice> tmp_lattices;
    std::map<unsigned int, Pin> tmp_pins;

    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    // Load XML file and put its contents in property tree.
    // No namespace qualification is needed, because of Koenig
    // lookup on the second argument. If reading fails, exception
    // is thrown.
    read_xml (filename, pt, boost::property_tree::xml_parser::trim_whitespace);
    std::string bin; // bin is a container to save temporal data
    // Dimension --------------------------------------------------------------
    tmp_dim = pt.get<unsigned int> ("geometry.<xmlattr>.dim", 0);
    Assert(tmp_dim == 1 or tmp_dim == 2 or tmp_dim == 3,
        dealii::ExcMessage("dimension should be 1, 2 or 3."));

    // core --------------------------------------------------------------
    {
    std::string name = pt.get<std::string> ("geometry.core.<xmlattr>.name", "Core");
    //unsigned int n_ref = pt.get<unsigned int> ("geometry.core.n_ref", 0);

    std::vector<unsigned int> n_nodes;
    bin = pt.get<std::string> ("geometry.core.nnodes");
    Forest::string_to_vector (bin, n_nodes);
    bin.clear();

    std::vector<std::vector<double> > length;
    bin = pt.get<std::string> ("geometry.core.length");
    Forest::string_to_vector (bin, length);
    bin.clear();

    std::vector<std::vector<unsigned int> > compositions;
    bin = pt.get<std::string> ("geometry.core.components");
    Forest::string_to_vector (bin, compositions);
    bin.clear();

    std::vector<unsigned int> bcs;
    bin = pt.get<std::string> ("geometry.core.boundary");
    Forest::string_to_vector (bin, bcs);
    bin.clear();

    tmp_core.set (n_nodes, length, compositions, bcs, name);
    }

    set(tmp_dim, tmp_core, tmp_lattices, tmp_pins);

    // assemblies --------------------------------------------------------
    for (const std::pair<std::string, ptree> & v : pt.get_child (
        "geometry.lattices"))
    {
      if (v.first == "lattice")
      {
        std::string name = v.second.get<std::string> ("<xmlattr>.name",
            "Lattice");
        unsigned int id = v.second.get<unsigned int> ("<xmlattr>.id");

        std::vector<unsigned int> n_nodes;
        bin = v.second.get<std::string> ("nnodes");
        Forest::string_to_vector (bin, n_nodes);
        bin.clear();

        std::vector<std::vector<unsigned int> > components;
        bin = v.second.get<std::string> ("components");
        Forest::string_to_vector (bin, components);
        bin.clear();

        std::vector<std::vector<double> > length;
        bin = v.second.get<std::string> ("length");
        Forest::string_to_vector (bin, length);
        bin.clear ();

        std::vector<double> water_gap;
        bin = v.second.get<std::string> ("water_gap", "0. 0. 0. 0.;");
        Forest::string_to_vector (bin, water_gap);
        bin.clear ();

        std::string latt_type = v.second.get<std::string> ("<xmlattr>.type");

        Lattice lattice;
        lattice.set (id, n_nodes, components,length, water_gap, latt_type, name);


        m_lattices.insert (
            std::pair<unsigned int, Lattice> (lattice.m_id, lattice));
      }
    }

    // pins --------------------------------------------------------------
    //for (auto v : pt.get_child("geometry.pins"))
    for (const std::pair<std::string, ptree> & v : pt.get_child (
        "geometry.pins"))
    {
      if (v.first == "pin")
      {
        Pin pin;

        unsigned int id = v.second.get<unsigned int> ("<xmlattr>.id");
        std::string pin_type = v.second.get<std::string> ("<xmlattr>.type");

        std::vector<unsigned int> materials;
        bin = v.second.get<std::string> ("materials");
        Forest::string_to_vector (bin, materials);
        bin.clear ();

        double fuel_radius = v.second.get<double> ("fuel_radius", 0.0);
        std::string name = v.second.get<std::string> ("<xmlattr>.name", "Pin");

        pin.set (id, materials, pin_type, fuel_radius , name);

        m_pins.insert (std::pair<unsigned int, Pin> (pin.m_id, pin));
      }
    }

    //set(tmp_dim, tmp_core, tmp_lattices, tmp_pins);
  }

  void
  InputGeom::save_xml (const std::string &filename)
  {
    std::ofstream out(filename);
    out << show();
    out.close();
  }

  std::string
  InputGeom::show ()
  {
    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    std::string bin; // bin is a container to save temporal data
    //std::string indent ("\t");
    std::string indent ("");

    pt.put ("geometry.<xmlattr>.dim", m_dim);
    // core --------------------------------------------------------------
    pt.add ("geometry.core", "");
    pt.put ("geometry.core.<xmlattr>.name", m_core.m_name);

    Forest::vector_to_string (m_core.m_n_nodes, bin);
    pt.put ("geometry.core.nnodes", bin);
    bin.clear ();

    Forest::vector_to_string (m_core.m_length, bin, indent, 3);
    pt.put ("geometry.core.length", bin);
    bin.clear ();

    Forest::vector_to_string (m_core.m_components, bin, indent, 3);
    pt.put ("geometry.core.components", bin);
    bin.clear ();

    Forest::vector_to_string (m_core.m_bcs, bin, indent, 3);
    pt.put ("geometry.core.boundary", bin);
    bin.clear ();

    // assemblies --------------------------------------------------------
    for (typename std::map<unsigned int, Lattice>::iterator lattice_it =
        m_lattices.begin (); lattice_it != m_lattices.end (); ++lattice_it)
    {
      // here we generate the node
      ptree & node = pt.add ("geometry.lattices.lattice", "");

      node.put ("<xmlattr>.id", lattice_it->first);
      Lattice &lattice = lattice_it->second;
      node.put ("<xmlattr>.name", lattice.m_name);

      Forest::vector_to_string (lattice.m_n_nodes, bin);
      node.put ("nnodes", bin);
      bin.clear ();

      Forest::vector_to_string (lattice.m_components, bin, indent, 4);
      node.put ("components", bin);
      bin.clear ();

      Forest::vector_to_string (lattice.m_length, bin, indent, 4);
      node.put ("length", bin);
      bin.clear ();

      if (lattice.m_latt_type == LattType::GRID)
        bin = "GRID";
      else if (lattice.m_latt_type == LattType::PINMAP)
        bin = "PINMAP";
      node.put ("type", bin);
      bin.clear();

      Forest::vector_to_string (lattice.m_water_gap, bin);
      node.put ("water_gap", bin);
      bin.clear ();

    }

    // pins --------------------------------------------------------------
    for (typename std::map<unsigned int, Pin>::iterator pin_it =
        m_pins.begin (); pin_it != m_pins.end (); ++pin_it)
    {
      ptree & node = pt.add ("geometry.pins.pin", "");

      Pin &pin = pin_it->second;
      node.put ("<xmlattr>.id", pin.m_id);
      node.put ("<xmlattr>.name", pin.m_name);

      Forest::vector_to_string (pin.m_materials, bin);
      node.put ("materials", bin);
      bin.clear ();

      if (pin.m_pin_type == PinType::PIN or pin.m_pin_type == PinType::PINNEW
          or pin.m_pin_type == PinType::PINBOX)
      {
        node.put ("fuel_radius", pin.m_fuel_radius);
      }
    }
    std::stringstream ss;
    write_xml (ss, pt, xml_writer_settings (' ', 1));
    return ss.str();
  }


  InputGeom
  InputGeom::generate_from_lattice (const unsigned int ass_id)
  {
    Core tmp_core;
    std::map<unsigned int, Lattice> tmp_lattices = {{ass_id, m_lattices[ass_id]}};
    //tmp_lattices.insert (std::pair<unsigned int, Lattice> (ass_id, m_lattices[ass_id]));

    InputGeom newgeo;
    newgeo.set(m_dim, tmp_core, tmp_lattices, m_pins);
    return newgeo;
  }

  void
  InputGeom::check_core () const
  {
    // Dimensions should be consistent

    unsigned int dim = m_dim;

    // more dimensions are not ready yet
    Assert(dim == 1 or dim == 2 or dim == 3,
        dealii::ExcMessage("dimension should be 1 or 2"));

    // First we check the dimensions with respect the number of nodes
    for (unsigned int i = dim; i < m_core.m_n_nodes.size (); ++i)
      Assert(1 == m_core.m_n_nodes[i],
          dealii::ExcMessage("n_nodes for dimensions not used should be one"));

    // Now we chek that the number of nodes per direction agree with the lengths
    for (unsigned int i = 0; i < dim; ++i)
      Assert(m_core.m_length[i].size () == m_core.m_n_nodes[i],
          dealii::ExcMessage("lengths vs n_nodes per axis error"));

    // Now we chek the number of nodes per axis versus the components
    if (dim == 1)
    {
      Assert(m_core.m_indices[0][0].size () == m_core.m_n_nodes[0],
          dealii::ExcMessage("lengths vs n_nodes per axis error"));
    }
    else
    {
      Assert(m_core.m_indices[0][0].size () == m_core.m_n_nodes[0],
          dealii::ExcMessage("lengths vs n_nodes per axis error"));
      for (unsigned int i = 0; i < m_core.m_indices[0].size (); ++i)
        Assert(m_core.m_indices[0][i].size () == m_core.m_n_nodes[dim - 2],
            dealii::ExcMessage("lengths vs n_nodes per axis error"));
    }

    // Now the boundary conditions
    Assert(m_core.m_bcs.size () == dim * 2,
        dealii::ExcMessage("boundary vs dimension*2 error"));

  } //end InputGeom::check_core()

  void
  InputGeom::check_lattices () const
  {
    unsigned int dim = m_dim;

    for (typename std::map<unsigned int, Lattice>::const_iterator lattice_it =
        m_lattices.begin (); lattice_it != m_lattices.end (); ++lattice_it)
    {
      const Lattice &lattice = lattice_it->second;

      // First we check the dimensions with respect the number of nodes
      for (unsigned int i = dim; i < lattice.m_n_nodes.size (); ++i)
        Assert(1 == lattice.m_n_nodes[i],
            dealii::ExcMessage("n_nodes for dimensions not used should be one"));

      // Depending on the type of lattice, the number of elements or the water gap
      // should be consistent
      if (lattice.m_latt_type == LattType::GRID)
      {
        for (unsigned int i = 0; i < dim; ++i)
          Assert(lattice.m_length[i].size () == lattice.m_n_nodes[i],
              dealii::ExcMessage("lengths vs n_nodes per axis error"));
      }
      else if (lattice.m_latt_type == LattType::PINMAP)
      {
        /*Assert(lattice->water_gap.size () == 2 * dim,
         dealii::ExcMessage("water_gap.size() vs dimension*2 error"));*/
      }
      else
      {
        Assert(false, dealii::ExcMessage("Latt_type not recognised"));
      }

      // Now we chek the number of nodes per axis versus the components
      if (dim == 1)
      {
        Assert(lattice.m_components.size () == dim,
            dealii::ExcMessage("number of lines in components is not valid for 1D"));
        Assert(lattice.m_components[0].size () == lattice.m_n_nodes[dim - 1],
            dealii::ExcMessage("components vs n_nodes per axis error"));
      }
      else
      {
        Assert(lattice.m_components.size () == lattice.m_n_nodes[dim - 1],
            dealii::ExcMessage("components vs n_nodes per axis error"));

        for (unsigned int i = 0; i < lattice.m_components.size (); ++i)
        {
          Assert(lattice.m_components[i].size () == lattice.m_n_nodes[dim - 2],
              dealii::ExcMessage("lattice: components vs n_nodes in axis " + Forest::num_to_string ( i) + " error\n"));
        }
      }
    }
  } // end InputGeom::check_lattices()

  void
  InputGeom::check () const
  {
    //save ("xml_test");

    this->check_core ();
    this->check_lattices ();
  } //end InputGeom::check()

} // end of namespace namespace Forest

