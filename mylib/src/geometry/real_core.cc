#include "geometry/real_core.h"

#include "input/input.h"
#include "input/input_geom.h"
#include "geometry/path_to_gridgenerator_dealii.h"
#include "utils/forest_utils_dealii.h"
#include "geometry/real_pin.h"
#include "geometry/real_assembly.h"

#include "deal.II/grid/tria.h"
#include "deal.II/base/point.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

namespace Forest
{
  using namespace dealii;

  template <int dim>
  RealCore<dim>::RealCore (const Input & input_data,
                           const double correction_factor)
      : m_input_data (input_data),
        m_core (input_data.mp_geom.m_core),
        m_origin (Point<dim> ()),
        m_ready_to_build_mesh (false),
        m_correction_factor (correction_factor),
        m_n_assemblies (0),
        m_n_pins (0)
  {
    this->build_data ();
  }

  template <int dim>
  void
  RealCore<dim>::build_data ()
  {
    this->build_assemblies ();
    this->build_additional_data ();
  }

  template <int dim>
  void
  RealCore<dim>::build_assemblies ()
  {
    // some cleaning to play safe
    m_real_pins.clear ();

    // prepare vector with the y component of the origin
    std::vector<double> z_origin (m_core.m_n_nodes[2]);
    z_origin[0] = 0.0;
    for (unsigned int k = 1; k < m_core.m_n_nodes[2]; ++k)
    {
      z_origin[k] = z_origin[k - 1] + m_core.m_length[2][k - 1];
    }

    // prepare vector with the y component of the origin
    std::vector<double> y_origin (m_core.m_n_nodes[1]);
    y_origin[0] = 0.0;
    for (unsigned int j = 1; j < m_core.m_n_nodes[1]; ++j)
    {
      y_origin[j] = y_origin[j - 1] + m_core.m_length[1][j - 1];
    }

    // prepare vector with the x component of the origin
    std::vector<double> x_origin (m_core.m_n_nodes[0]);
    x_origin[0] = 0.0;
    for (unsigned int i = 1; i < m_core.m_n_nodes[0]; ++i)
    {
      x_origin[i] = x_origin[i - 1] + m_core.m_length[0][i - 1];
    }

    // fill the vector of assemblies with the data of each assembly
    for (unsigned int j = 0; j < m_core.m_n_nodes[1]; ++j)
    {
      for (unsigned int i = 0; i < m_core.m_n_nodes[0]; ++i)
      {
        if (m_core.m_indices[0][j][i] != m_empty_assembly)
        {
          RealAssembly<dim> real_assembly;

          unsigned int ind = m_core.m_indices[0][j][i];
          unsigned int comp_ind = m_core.m_components[ind];
          const Lattice & lattice = m_input_data.mp_geom.m_lattices.at(comp_ind);
          Point<dim> latt_origin (
              new_point<dim> (x_origin[i], y_origin[j], 0.0));
          real_assembly.init (lattice, latt_origin, m_correction_factor);

          // this is to plot the information of the particular assembly
          //real_assembly.info ();

          const std::map<unsigned int, Pin> & pins = m_input_data.mp_geom.m_pins;
          real_assembly.build_data (pins, m_real_pins);

          m_real_assemblies.push_back (real_assembly);

          /*
           unsigned int non_empty_pins = data_pins.size ();
           std::vector<unsigned int> pins_this_latt (non_empty_pins,
           map_latt_index);
           pins_to_assembly.insert (pins_to_assembly.end (),
           pins_this_latt.begin (), pins_this_latt.end ());
           */
        } //if
      }
    }

    m_ready_to_build_mesh = true;
  }

  template <int dim>
  void
  RealCore<dim>::build_additional_data ()
  {
    m_n_assemblies = m_real_assemblies.size ();

    m_n_pins = 0;
    for (unsigned int i = 0; i < m_real_assemblies.size (); i++)
    {
      m_n_pins += m_real_assemblies[i].get_n_pins ();
    }

    pins_per_assembly.clear ();
    for (unsigned int i = 0; i < m_real_assemblies.size (); i++)
    {
      pins_per_assembly.push_back (m_real_assemblies[i].get_n_pins ());
    }

    pins_to_assembly.clear ();
    for (unsigned int i = 0; i < m_real_assemblies.size (); i++)
    {
      for (unsigned int j = 0; j < m_real_assemblies[i].get_n_pins (); j++)
      {
        pins_to_assembly.push_back (i);
      }
    }
  }

  template <int dim>
  void
  RealCore<dim>::build_mesh (Triangulation<dim, dim> & tria)
  {

    Assert(m_n_assemblies == m_real_assemblies.size (),
        ExcMessage("Wrong number of real assemblies."));

    std::vector<Triangulation<dim, dim> > m_tria_latt (m_n_assemblies);
    for (unsigned int i = 0; i < m_real_assemblies.size (); ++i)
    {
      m_real_assemblies[i].build_mesh (m_real_pins, m_tria_latt[i]);
    }

    GridGenerator::merge_triangulations (m_tria_latt, tria);
  }

  template class RealCore<1> ;
  template class RealCore<2> ;
  template class RealCore<3> ;

} // end of namespace Forest
