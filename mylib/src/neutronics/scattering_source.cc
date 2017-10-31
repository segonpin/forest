/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_source.cc
 * @brief  Implementation of class template ScatteringSource
 */

#include "neutronics/scattering_source.h"

#include "algebra/sn_vector.h"        // for SnVector

#include "deal.II/lac/block_vector.h" // for BlockVector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  void
  ScatteringSource<dim>::down_scatt_add (SnVector & source,
                                         const SnVector & phi,
                                         const unsigned int g)
  {
    for (unsigned int h = 0; h < g; ++h)
    {
      vmult_add_g_h (source.m_data[g], phi.m_data[h], g, h);
    }
  }

  template <int dim>
  void
  ScatteringSource<dim>::up_scatt_add (SnVector & source,
                                       const SnVector & phi,
                                       const unsigned int g)
  {
    for (unsigned int h = g + 1; h < phi.m_data.size (); ++h)
    {
      vmult_add_g_h (source.m_data[g], phi.m_data[h], g, h);
    }
  }

  template <int dim>
  void
  ScatteringSource<dim>::in_scatt (SnVector & source,
                                   const SnVector & phi,
                                   const unsigned int g)
  {
    in_scatt (source.m_data[g], phi.m_data[g], g);
  }

  template <int dim>
  void
  ScatteringSource<dim>::in_scatt (BlockVector<double> & source,
                                   const BlockVector<double> & phi,
                                   const unsigned int g)
  {
    source = 0;
    vmult_add_g_h (source, phi, g, g);
  }

  template <int dim>
  void
  ScatteringSource<dim>::in_scatt_add (SnVector & source,
                                       const SnVector & phi,
                                       const unsigned int g)
  {
    vmult_add_g_h (source.m_data[g], phi.m_data[g], g, g);
  }

  template <int dim>
  void
  ScatteringSource<dim>::vmult_add_g_h (SnVector & source,
                                        const SnVector & phi,
                                        const unsigned int g,
                                        const unsigned int h)
  {
    vmult_add_g_h (source.m_data[g], phi.m_data[h], g, h);
  }

  template class ScatteringSource<1> ;
  template class ScatteringSource<2> ;
  template class ScatteringSource<3> ;

} // end of namespace Forest

