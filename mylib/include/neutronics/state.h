/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/state.h
 * @brief  State class template declarations
 */

#ifndef FOREST_STATE_H
#define FOREST_STATE_H

#include "input/input_fwd.h"            // For class Input&
#include "algebra/sn_vector.h"          // for SnVector
#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&

#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/fe/fe_dgq.h"          // for FE_DGQ
#include "deal.II/fe/mapping_q1.h"      // for MappingQ1

namespace Forest
{
  using namespace dealii;

  /**
   * @brief State vectors for the system
   * @ingroup ForestNeutronics
   */
  template <int dim>
  class State
  {
  public:
    /**
     * @brief Constructor.
     * @param data
     * @param geom
     */
    State (const Input &data,
           const ProbGeom<dim> &geom);

    /** @brief Indicating the level of verbose: 0 = none; 1 = everyting. */
    bool m_verbose;

    /** @brief Input data */
    const Input &mp_data;
    /** @brief Geometry */
    const ProbGeom<dim> &mp_geom;

  private:
    /** @todo Why do we need this? */
    const MappingQ1<dim> m_mapping;

    FE_DGQ<dim> m_fe;

    DoFHandler<dim, dim> m_dof_handler;
    unsigned int m_n_groups;
    unsigned int m_n_leg_mom; /** @todo move this to the materials input? */
    unsigned int m_n_elem;
    double m_keff;

  public:

    /** @brief Scalar flux */
    SnVector m_scalar_flux;

    /** @brief Container for the angular flux. */
    SnVector m_psi;

    /**
     * @brief Return the memory consumption of an object of this class.
     * @return The memory consumption
     * @todo Add mapping.memory_consumption() and fe.memory_consumption()?
     */
    double
    memory_consumption () const;

    /** @todo document me */
    unsigned int
    get_n_groups () const
    {
      return m_n_groups;
    }

    /** @todo document me */
    unsigned int
    get_n_elem () const
    {
      return m_n_elem;
    }

    /** @todo document me */
    unsigned int
    get_n_leg_mom () const
    {
      return m_n_leg_mom;
    }

    /** @todo document me */
    const FiniteElement<dim, dim> &
    get_fe () const
    {
      return m_fe;
    }

    /** @todo document me */
    const MappingQ1<dim> &
    get_mapping () const
    {
      return m_mapping;
    }

    /** @todo document me */
    const DoFHandler<dim, dim> &
    get_dof_handler () const
    {
      return m_dof_handler;
    }

    /**
     * @brief Set the k-effective with a given value.
     * @param keff
     */
    void
    set_keff (double keff)
    {
      m_keff = keff;
    }

    /**
     * @brief Return the keff associated with this state.
     * @return keff
     */
    double
    get_keff () const
    {
      return m_keff;
    }

  };
} // end of namespace Forest

#endif /* FOREST_STATE_H */
