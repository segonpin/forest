/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/manager.h
 * @brief  State class template declarations
 */

#ifndef FOREST_MANAGER_H
#define FOREST_MANAGER_H

#include "input/input.h"            // For class Input&
#include "algebra/sn_vector.h"          // for SnVector
#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&

namespace Forest
{
  using namespace dealii;

  /**
   * @brief State vectors for the system
   * @ingroup ForestNeutronics
   */
  template <int dim>
  class Manager
  {
  public:
    /**
     * @brief Constructor.
     * @param data
     * @param geom
     */
    Manager (Input & data)
    : mp_data(data)
    {
    };

    void
    set_lattice_problems()
    {
      std::string path = mp_data.mp_settings.path;
      std::string filename_no_ext = mp_data.mp_settings.problem_name;

      // we have just one core
      core_data = mp_data;

      // we have many lattices
      unsigned int n_lattices = mp_data.mp_geom.m_lattices.size();
      ass_data.resize(n_lattices);
      for (unsigned int i = 0; i < n_lattices; ++i)
      {
        ass_data[i].mp_settings = mp_data.mp_settings;
        ass_data[i].mp_mat = mp_data.mp_mat;
        ass_data[i].mp_geom = mp_data.mp_geom.generate_from_lattice(i);
      }
    }

    /** @brief Input data */
    Input & mp_data;

    Input core_data;
    std::vector<Input> ass_data;

    /**
     * @brief Return the memory consumption of an object of this class.
     * @return The memory consumption
     * @todo Add mapping.memory_consumption() and fe.memory_consumption()?
     */
    double
    memory_consumption () const
    {
      double memory_consumption = 0.0;
      return memory_consumption;
    };


  };
} // end of namespace Forest

#endif /* FOREST_MANAGER_H */
