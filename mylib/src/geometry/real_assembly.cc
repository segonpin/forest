/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   geometry/real_assembly.cc
 * @brief  Implementation of RealAssembly functions
 */
#include "geometry/real_assembly.h"

#include "geometry/path_to_gridgenerator_dealii.h"
#include "geometry/real_pin.h"
#include "utils/forest_utils_dealii.h"

#include "deal.II/grid/tria.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <iostream>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  RealAssembly<dim>::RealAssembly ()
      : m_n_nodes (),
        m_components (),
        m_type (),
        m_length (),
        m_origin (),
        m_pins_ini (0),
        m_pins_end (0),
        m_n_pins (0),
        m_correction_factor (1.0),
        m_initialized (false),
        m_ready_to_build_mesh (false)
  {
  }

  template <int dim>
  void
  RealAssembly<dim>::init (const Forest::Lattice & lattice,
                           Point<dim> origin,
                           const double correction_factor)
  {
    this->set_lattice_data (lattice);
    this->set_origin (origin);
    this->set_correction_factor (correction_factor);
    m_initialized = true;
  }

  template <int dim>
  void
  RealAssembly<dim>::build_data (const std::map<unsigned int, Pin> & pins,
                                 std::vector<RealPin<dim> > & real_pins)
  {
    this->build_pins (pins, real_pins);
    this->build_additional_data ();
  }

  template <int dim>
  void
  RealAssembly<dim>::set_lattice_data (const Forest::Lattice & lattice)
  {
    m_n_nodes = lattice.m_n_nodes;
    m_components = lattice.m_components;
    m_type = lattice.m_latt_type;
    /*if (m_type == LattType::PINMAP)
      m_pitch = new_point<dim> (lattice.m_pitch[0], lattice.m_pitch[0],
          lattice.m_pitch[0]);
    else*/
    m_length = lattice.m_length;
  }

  template <int dim>
  void
  RealAssembly<dim>::build_pins (const std::map<unsigned int, Pin> & pins,
                                 std::vector<RealPin<dim> > &real_pins)
  {
    Assert(m_initialized, ExcMessage("Object not initialized yet."));
    m_pins_ini = real_pins.size ();
    if (m_type == Forest::LattType::PINMAP)
      set_assembly_pin_map (pins, real_pins);
    else if (m_type == Forest::LattType::GRID)
      set_assembly_grid (real_pins);
    else
      Assert(false, ExcNotImplemented());
    m_pins_end = real_pins.size ();
    m_ready_to_build_mesh = true;
  }

  template <int dim>
  void
  RealAssembly<dim>::build_additional_data ()
  {
    m_n_pins = m_pins_end - m_pins_ini;
  }

  template <int dim>
  void
  RealAssembly<dim>::set_assembly_pin_map (const std::map<unsigned int, Pin> & pins,
                                           std::vector<RealPin<dim> > &real_pins)
  {
    // prepare vector with the z component of the origin
    std::vector<double> z_origin (m_n_nodes[2]);
    z_origin[0] = 0.0;
    for (unsigned int k = 1; k < m_n_nodes[2]; ++k)
      z_origin[k] = z_origin[k - 1] + m_length[2][k - 1];

    // prepare vector with the y component of the origin
    std::vector<double> y_origin (m_n_nodes[1]);
    y_origin[0] = 0.0;
    for (unsigned int j = 1; j < m_n_nodes[1]; ++j)
      y_origin[j] = y_origin[j - 1] + m_length[1][j - 1];

    // prepare vector with the x component of the origin
    std::vector<double> x_origin (m_n_nodes[0]);
    x_origin[0] = 0.0;
    for (unsigned int i = 1; i < m_n_nodes[0]; ++i)
      x_origin[i] = x_origin[i - 1] + m_length[0][i - 1];

    // we generate the information for the assemblies
    for (unsigned int j = 0; j < m_n_nodes[1]; ++j)
      for (unsigned int i = 0; i < m_n_nodes[0]; ++i)
        if (m_components[j][i] != m_empty_pin)
        {
          Point<dim> pin_origin = m_origin
              + Forest::new_point<dim> (x_origin[i], y_origin[j], 0.0);

          //RealPin<dim> real_pin;
          //const Forest::Pin & pin = pins.at(m_components[j][i]);
          RealPin<dim> real_pin;
          real_pin.init (pins.at(m_components[j][i]), m_correction_factor);
          if(dim==1)
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], 0., 0.));
          }
          else if(dim==2)
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], m_length[1][j], 0.));
          }
          else
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], m_length[1][j], m_length[2][j]));
          }
          real_pin.set_center (pin_origin);
          // real_pin.info ();
          real_pins.push_back (real_pin);
        }

    /*
    for (unsigned int j = 0; j < m_n_nodes[1]; ++j)
      for (unsigned int i = 0; i < m_n_nodes[0]; ++i)
        if (m_components[j][i] != m_empty_pin)
        {
          Point<dim> pin_origin = m_origin
              + new_point<dim> (i * m_pitch[0], j * m_pitch[0], 0.0);
          const Forest::Pin & pin = pins.at(m_components[j][i]);
          RealPin<dim> real_pin;
          real_pin.init (pin, m_correction_factor);
          Point<dim> pitch_aux (
              new_point<dim> (m_pitch[0], m_pitch[0], m_pitch[0]));
          real_pin.set_pitch (pitch_aux);
          real_pin.set_center (pin_origin);
          //real_pin.info ();
          real_pins.push_back (real_pin);
        }
   */
  }

  template <int dim>
  void
  RealAssembly<dim>::set_assembly_grid (std::vector<RealPin<dim> > &real_pins)
  {
    // prepare vector with the z component of the origin
    std::vector<double> z_origin (m_n_nodes[2]);
    z_origin[0] = 0.0;
    for (unsigned int k = 1; k < m_n_nodes[2]; ++k)
      z_origin[k] = z_origin[k - 1] + m_length[2][k - 1];

    // prepare vector with the y component of the origin
    std::vector<double> y_origin (m_n_nodes[1]);
    y_origin[0] = 0.0;
    for (unsigned int j = 1; j < m_n_nodes[1]; ++j)
      y_origin[j] = y_origin[j - 1] + m_length[1][j - 1];

    // prepare vector with the x component of the origin
    std::vector<double> x_origin (m_n_nodes[0]);
    x_origin[0] = 0.0;
    for (unsigned int i = 1; i < m_n_nodes[0]; ++i)
      x_origin[i] = x_origin[i - 1] + m_length[0][i - 1];

    // we generate the information for the assemblies
    for (unsigned int j = 0; j < m_n_nodes[1]; ++j)
      for (unsigned int i = 0; i < m_n_nodes[0]; ++i)
        if (m_components[j][i] != m_empty_pin)
        {
          Point<dim> pin_origin = m_origin
              + Forest::new_point<dim> (x_origin[i], y_origin[j], 0.0);

          RealPin<dim> real_pin;
          real_pin.init (m_components[j][i]);
          if(dim==1)
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], 0., 0.));
          }
          else if (dim ==2)
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], m_length[1][j], 0.));
          }
          else
          {
            real_pin.set_pitch (
                new_point<dim> (m_length[0][i], m_length[1][j], m_length[2][j]));
          }
          real_pin.set_center (pin_origin);
          // real_pin.info ();
          real_pins.push_back (real_pin);
        }
  }

  template <int dim>
  void
  RealAssembly<dim>::build_mesh (std::vector<RealPin<dim> > &real_pins,
                                 Triangulation<dim, dim> & tria)
  {
    std::vector<Triangulation<dim, dim> > m_tria_pins (m_pins_end - m_pins_ini);
    for (unsigned int i = m_pins_ini; i < m_pins_end; ++i)
      real_pins[i].build_mesh (m_tria_pins[i - m_pins_ini]);

    // We merge the pins in a new assembly
    GridGenerator::merge_triangulations (m_tria_pins, tria);
  }

  /**
   * @todo I have to complete this
   */
  template <int dim>
  void
  RealAssembly<dim>::info () const
  {
    std::cout << "Information for assembly of type:";
    if (m_type == Forest::LattType::PINMAP)
      std::cout << "LattType::PINMAP" << std::endl;
    else if (m_type == Forest::LattType::GRID)
      std::cout << "LattType::GRID" << std::endl;
    else
      Assert(false, ExcMessage("Wrong type of assembly."));

    std::cout << "m_origin = " << m_origin << std::endl;

    std::cout << std::endl;
  }

  template class RealAssembly<1> ;
  template class RealAssembly<2> ;
  template class RealAssembly<3> ;

} // end of namespace Forest
