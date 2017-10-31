#include "geometry/real_pin.h"

#include "geometry/path_to_gridgenerator_dealii.h"
#include "input/input_geom.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage
#include "deal.II/base/point.h"

#include <iostream>
#include <algorithm>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  RealPin<dim>::RealPin ()
      : m_pitch (),
        m_center (),
        m_type (Forest::PinType::BOX),
        m_materials (),
        m_radius (0.0)
  {
  }

  template <int dim>
  void
  RealPin<dim>::init (const Forest::Pin & pin,
                      const double correction_factor)
  {
    m_type = pin.m_pin_type;
    m_materials = pin.m_materials;
    if (m_type == Forest::PinType::PIN or m_type == Forest::PinType::PINNEW)
    {
      m_radius = pin.m_fuel_radius * correction_factor;
    }
    else if (m_type == Forest::PinType::PINBOX)
    {
      // this is to recover the fact that when making the pin, we try to
      // change the radius to match the circle
      m_radius = pin.m_fuel_radius * 2./std::sqrt (2.0);
    }
  }

  template <int dim>
  void
  RealPin<dim>::init (const unsigned int material_id)
  {
    m_type = Forest::PinType::BOX;
    m_materials.clear ();
    m_materials.push_back (material_id);
  }

  template <int dim>
  bool
  RealPin<dim>::is_inside (const Point<dim> p) const
  {
    bool inside = true;
    for (unsigned int d = 0; d < dim; ++d)
    {
      if (std::abs (p (d) - m_center (d)) >= m_pitch (d) / 2.0)
      {
        inside = false;
      }
    }
    return inside;
  }

  template <int dim>
  double
  RealPin<dim>::norm_inf_xy (const Point<dim> & point) const
  {
    double norm = 0;
    for (unsigned int i = 0; i < std::min (2, dim); ++i)
    {
      norm = std::max (norm, std::abs (m_center (i) - point (i)));
    }
    return norm;
  }

  template <int dim>
  void
  RealPin<dim>::build_mesh (Triangulation<dim, dim> &tria)
  {
    if (m_type == Forest::PinType::PIN)
    {
      GridGenerator::pin_cell (tria, m_center, m_pitch, m_radius, m_materials);
    }
    else if (m_type == Forest::PinType::PINNEW)
    {
      GridGenerator::pin_cell2 (tria, m_center, m_pitch, m_radius, m_materials);
    }
    else if (m_type == Forest::PinType::PINBOX)
    {
      GridGenerator::pin_cell (tria, m_center, m_pitch, m_radius, m_materials);
    }
    else if (m_type == Forest::PinType::BOX)
    {
      std::cout << " m_pitch = " << m_pitch << std::endl;
      GridGenerator::pin_box (tria, m_center, m_pitch, m_materials);
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }
  }

  template <int dim>
  void
  RealPin<dim>::info () const
  {

    std::cout << "Information for pin of type:";
    if (m_type == Forest::PinType::PIN)
    {
      std::cout << "PinType::PIN" << std::endl;
    }
    else if (m_type == Forest::PinType::PINNEW)
    {
      std::cout << "PinType::PINNEW" << std::endl;
    }
    else if (m_type == Forest::PinType::PINBOX)
    {
      std::cout << "PinType::PINBOX" << std::endl;
    }
    else if (m_type == Forest::PinType::BOX)
    {
      std::cout << "PinType::BOX" << std::endl;
    }
    else
    {
      Assert(false, ExcMessage("Wrong pin type."));
    }

    std::cout << "m_center = " << m_center << std::endl;

    std::cout << "m_pitch = " << m_pitch << std::endl;
    std::cout << "m_materials =";
    for (unsigned int i = 0; i < m_materials.size (); ++i)
    {
      std::cout << " " << m_materials[i];
    }
    std::cout << std::endl;
    if (m_type == Forest::PinType::PIN or m_type == Forest::PinType::PINBOX)
    {
      std::cout << "m_radius = " << m_radius << std::endl;
    }
    std::cout << std::endl;
  }

  template class RealPin<1> ;
  template class RealPin<2> ;
  template class RealPin<3> ;

} // end of namespace Forest
