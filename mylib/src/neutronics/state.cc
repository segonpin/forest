/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/state.cc
 * @brief  Implementation of class template State
 */

#include "neutronics/state.h"

#include "geometry/prob_geom.h"
#include "input/input.h"
#include "utils/forest_utils_logstream.h"
#include "utils/forest_utils_memory.h"
#include "utils/forest_utils_timer.h"

#include "deal.II/base/timer.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/grid/tria.h"

namespace Forest
{
  template <int dim>
  State<dim>::State (const Input &data,
                     const ProbGeom<dim> &geom)
      : m_verbose (false),
        mp_data (data),
        mp_geom (geom),
        m_mapping (),
        m_fe (data.mp_settings.get_fe_degree ()),
        m_dof_handler (geom.m_tria)
  {
    if (m_verbose)
      log_print_text ("Generating the state object", 1);

    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart ();

    /* We initialize some data. */
    m_n_groups = mp_data.mp_mat.get_n_groups ();
    m_n_leg_mom = 1;
    m_n_elem = m_fe.dofs_per_cell * mp_geom.m_tria.n_active_cells ();
    m_keff = 1.0;

    /* Distribute dofs. */
    m_dof_handler.distribute_dofs (m_fe);

    // Stop the timer and report the time
    timer.stop ();
    StatsTimer::instance ().add ("Generating State",
        timer.wall_time ());
    // Report the memory
    StatsMemory::instance ().add ("State", this->memory_consumption ());
  }

  template <int dim>
  double
  State<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption += m_dof_handler.memory_consumption ();
    memory_consumption += m_scalar_flux.memory_consumption ();
    memory_consumption += m_psi.memory_consumption ();
    return memory_consumption;
  }

  template class State<1> ;
  template class State<2> ;
  template class State<3> ;

} // end of namespace Forest
