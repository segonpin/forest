/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/solver_eq.cc
 * @brief  Implementation of class templates SolverEq and EqSolverWG
 */

#include "algebra/solver_eq.h"

#include "algebra/equation.h"
#include "algebra/sn_vector.h"             // for SnVector
#include "neutronics/state.h"
#include "input/input.h"

#include "utils/forest_utils_logstream.h" // for log_print_text
#include "utils/forest_utils_memory.h"   // for StatsMemory
#include "utils/forest_utils_timer.h"    // for StatsTimer

#include "deal.II/base/timer.h"          // for Timer
#include "deal.II/lac/precondition.h"   // for PreconditionIdentity
#include "deal.II/lac/solver_control.h" // for SolverControl
#include "deal.II/lac/solver_gmres.h"   // for SolverGMRES, SolverFGMRES, etc
#include "deal.II/lac/solver_cg.h"   // for SolverCG
#include "deal.II/lac/solver_richardson.h"   // for SolverRICHARDSON
#include "deal.II/lac/solver_selector.h"
#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <vector>                       // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  SolverEq<dim>::SolverEq (State<dim> & state,
                           std::shared_ptr<Equation<dim> > & equation)
      : mp_state (state),
        m_vector_memory (8),
        m_max_it_i (mp_state.mp_data.mp_settings.get_inner_max_it ()),
        m_tol_i (mp_state.mp_data.mp_settings.get_inner_tol ())
  {
    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart (); // We restart the timer

    /* Set the control parameters. */
    std::string bin;
    bin = mp_state.mp_data.mp_settings.get_mg_solver();
    deallog << "mg_solver = " << mp_state.mp_data.mp_settings.get_mg_solver() << std::endl;
    m_use_gauss_seidel = bin.compare (std::string ("GS")) == 0;
    deallog << "m_use_gauss_seidel = " << m_use_gauss_seidel << std::endl;

    bin.clear(); bin = mp_state.mp_data.mp_settings.get_inner_solver();
    m_use_krylov = bin.compare (std::string ("krylov")) == 0;

    /* These constants are used to initialize some vectors. */
    m_n_leg_mom = mp_state.get_n_leg_mom ();
    m_n_elem = mp_state.get_n_elem ();
    m_n_groups = mp_state.get_n_groups ();

    /* Initialization of the counters for the sweeps. */
    m_total_sweeps = 0;
    m_total_sweeps_this_iteration = 0;

    /* Initialize auxiliary vectors. */
    m_sn_phi.reinit (m_n_elem);
    m_tmp_source.reinit (m_n_elem);
    m_aux_vec.reinit (m_n_leg_mom, m_n_elem);
    m_outer_source.reinit (m_n_leg_mom, m_n_elem);

    /* vector with the size of each block for flux and boundary conditions. */
    m_flux_bv_sizes = equation->get_flux_bv_sizes ();

    /*
    deallog << "m_flux_bv_sizes.size() = " << m_flux_bv_sizes.size() << std::endl;
    for (unsigned int i = 0; i < m_flux_bv_sizes.size(); ++i)
      deallog << "m_flux_bv_sizes[i] = " << m_flux_bv_sizes[i] << std::endl;
    */

    /* Declaring the rhs and the solution vectors. */
    m_rhs.reinit (m_flux_bv_sizes);
    m_rhs.collect_sizes ();

    m_mgrhs.reinit(m_n_groups, m_flux_bv_sizes);
    m_sol.reinit (m_flux_bv_sizes);
    m_sol.collect_sizes ();
    /* Auxiliary vector to store old solution to be used as next guess. */
    m_phi_old.reinit (m_n_groups, m_n_leg_mom, m_n_elem);
    /* Declaring auxiliary vectors. */
    m_sol_old.reinit (m_flux_bv_sizes);
    m_sol_old.collect_sizes ();

    // Stop the timer and report the time
    timer.stop ();
    StatsTimer::instance ().add ("SolverEq", timer.wall_time ());
    // Report the memory
    StatsMemory::instance ().add ("SolverEq", memory_consumption ());
  }

  template <int dim>
  void
  SolverEq<dim>::iterate (std::shared_ptr<Equation<dim> > & equation,
                          SnVector & phi_new,
                          const SnVector & fis_source_old,
                          const double /*outer_error*/) const
  {

    // Initialization of the counters for the sweeps.
    m_total_sweeps_this_iteration = 0;

    // auxiliary structure to efficiently implement the Gauss-Seidel
    std::vector<bool> skip_group (m_n_groups, false);

    // We prepare this for several iterations of the Gauss-Seidel, but the
    // proper way must be using Gauss-Seidel to precondition Krylov in a
    // multigroup problem
    //const unsigned int gs_iter = 1;

    unsigned int max_it = mp_state.mp_data.mp_settings.get_mg_max_it();
    double tol = mp_state.mp_data.mp_settings.get_mg_tol();;

    if(m_use_gauss_seidel)
    {
      for (unsigned int iter = 0;;iter++)
      {
        // phi_new contains the initial guess (actually the last solution)
        m_phi_old = phi_new;
        gauss_seidel_LDU (equation, phi_new,  m_phi_old, fis_source_old);

        m_phi_old.sadd(-1.,1.,phi_new);
        double correction_size = m_phi_old.l2_norm();
        if (m_phi_old.l2_norm() < tol || iter >= max_it)
        {
          deallog << "gs iters = " << iter << " and correction = " << correction_size << std::endl;
          break;
        }
      }
    }
    else
    {
      // First we prepare the multigroup right hand side
      for (unsigned int g = 0; g < m_n_groups; ++g)
      {
        //deallog << "preparing rhs for g = " << g << std::endl;
        //deallog << "m_mgrhs.m_data[g].n_blocks() = " << m_mgrhs.m_data[g].n_blocks() << std::endl;

        // Ther within group source starts with the fission using the old flux
        m_outer_source = fis_source_old.m_data[g];
        // We tell the energy group to the matrix vector product.
        equation->set_group (g);
        // First we prepare the right hand side
        equation->set_rhs (m_mgrhs.m_data[g], m_outer_source);
      }

      Multigroup<dim> multigroup_equation(equation);
      MGPrec<dim> multigroup_preconditioner(multigroup_equation);
      SolverControl solver_control (max_it, tol);
      SolverGMRES<SnVector > solver (solver_control, m_mg_vector_memory,
          SolverGMRES<SnVector >::AdditionalData (8));
      //solver.solve (multigroup_equation, phi_new, m_mgrhs, multigroup_preconditioner);
      solver.solve (multigroup_equation, phi_new, m_mgrhs, PreconditionIdentity ());
      deallog << "gmres iters = " << solver_control.last_step () << " and residual = " << solver_control.last_value () << std::endl;


    }

    m_total_sweeps += m_total_sweeps_this_iteration;
  }

  template <int dim>
  void
  SolverEq<dim>::gauss_seidel_LDU (std::shared_ptr<Equation<dim> > & equation,
      SnVector & phi_new,
      const SnVector & phi_old,
      const SnVector & fis_source_old) const
  {
      // Loop over the energy groups to solve the source iteration per group.
      for (unsigned int g = 0; g < m_n_groups; ++g)
      {
        /*if (skip_group[g])
        {
          continue;
        }*/
        // Ther within group source starts with the fission using the old flux
        m_outer_source = fis_source_old.m_data[g];
        // The old flux is the best guess for the new flux.
        phi_new.m_data[g] = phi_old.m_data[g];
        // Add up-scatt with the old flux to the external source.
        equation->up_add (m_outer_source, phi_old, g);
        // Add down-scatt. with new flux (if Gauss-Seidel) or old otherwise.
        equation->down_add (m_outer_source, phi_new, g);

        // We tell the energy group to the matrix vector product.
        equation->set_group (g);
        // First we prepare the right hand side
        equation->set_rhs (m_rhs, m_outer_source);
        // we add one sweep that is needed to generate the rhs
        ++m_total_sweeps_this_iteration;

        // Put the flux and the boundary conditions in one vector.
        m_sol = phi_old.m_data[g];

        // Use Krylov method or Richardson iteration.
        // const double new_tol = m_tol_i / outer_error;
        // const unsigned int new_max_it = m_max_it_i;
        /// @todo add new_tol, new_max_it?
        deallog << "g = " << g << "; ";

        // this represents the internal dofs, and we want the ones in the boundary


        equation->solve (m_sol, m_rhs);

        /*
        unsigned int n_elem = mp_state.get_n_elem();
        deallog << " m_rhs =";
        for (unsigned int i = n_elem ; i < m_rhs.size(); i++)
          deallog << " " << m_rhs[i];
        deallog  << std::endl;

        deallog << " m_sol =";
        for (unsigned int i = n_elem ; i < m_sol.size(); i++)
          deallog << " " << m_sol[i];
        deallog  << std::endl;
        deallog  << std::endl;
        */


        // Extract the solution and the boundary conditions from the one vector.
        phi_new.m_data[g] = m_sol;

        // Accumulate the within group iterations to the total iterations.
        m_total_sweeps_this_iteration += equation->get_tot_wg_it ();

        // If there is no upscattering for this groups and the previous ones,
        // we do not use them again in the Gauss Seidel.
        /*
        if (g == 0)
          if (equation->no_upscatt (g))
            skip_group[g] = true;
        else
          if (equation->no_upscatt (g) and skip_group[g - 1])
            skip_group[g] = true;
        */
      }
  }


  template <int dim>
  double
  SolverEq<dim>::get_sweeps_per_group () const
  {
    return double (m_total_sweeps_this_iteration) / double (m_n_groups);
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  template <int dim>
  EqSolverWG<dim>::EqSolverWG (State<dim> & state,
                               Equation<dim> & equation)
      : mp_state (state),
        mp_equation (equation),
        m_preconditioner (mp_equation), // Parsing the preconditioner
        m_vector_memory (8),
        m_max_it_i (mp_state.mp_data.mp_settings.get_inner_max_it ()),
        m_tol_i (mp_state.mp_data.mp_settings.get_inner_tol ())
  {
    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart (); // We restart the timer
    // Set the control parameters.
    m_use_krylov = true;
    m_use_richardson = false;
    // Initialization of the counters for the sweeps.
    m_total_wg_it = 0;
    // Stop the timer and report the time
    timer.stop ();
    StatsTimer::instance ().add ("EqSolverWG", timer.wall_time ());
    // Report the memory
    StatsMemory::instance ().add ("EqSolverWG", memory_consumption ());
  }

  template <int dim>
  void
  EqSolverWG<dim>::solve (BlockVector<double> & sol,
                          const BlockVector<double> & rhs,
                          const double tol_default,
                          const unsigned int max_it_default) const
  {
    const double tol = (tol_default < 0) ? m_tol_i : tol_default;
    const unsigned int max_it =
        (max_it_default == 0) ? m_max_it_i : max_it_default;

    m_total_wg_it = 0;
    // Use Krylov method or Richardson iteration.
    std::string solver_name = "gmres"; //"fgmres or gmres"
    if (m_use_krylov)
    {
      if (mp_equation.is_symmetric ())
        solver_name = "cg";
      else
        solver_name = "gmres";
    }
    else if (m_use_richardson)
    {
      // phi_new = [D*L_inv*M*S]*phi_old
      //         = [I-(I-D*L_inv*M*S)]*phi_old ;
      //         = [I-Eq.vmult]*phi_old ;
      solver_name = "richardson";
    }
    else
    {
      Assert(false, ExcMessage("Krylov or Richardson must be chosen."));
    }

    //--------------------------------------------------------------------------
    // Prepare gmres to solve the system.
    SolverControl solver_control (max_it, tol);
    //--------------------------------------------------------------------------
    // By default we use the FGMRES
    if (solver_name == "fgmres")
    {
      SolverFGMRES<BlockVector<double> > solver (solver_control,
          m_vector_memory,
          SolverFGMRES<BlockVector<double> >::AdditionalData (8));
      solver.solve (mp_equation, sol, rhs, m_preconditioner);
      deallog << "fgmres iters = " << solver_control.last_step () << std::endl;
    }
    //--------------------------------------------------------------------------
    // we use GMRES
    else if (solver_name == "gmres")
    {
      SolverGMRES<BlockVector<double> > solver (solver_control, m_vector_memory,
          SolverGMRES<BlockVector<double> >::AdditionalData (8));
      solver.solve (mp_equation, sol, rhs, m_preconditioner);
      deallog << "gmres iters = " << solver_control.last_step () << " and tol = " << solver_control.last_value () << std::endl;
    }
    //--------------------------------------------------------------------------
    // we use CG
    else if (solver_name == "cg")
    {
      SolverCG<BlockVector<double> > solver (solver_control, m_vector_memory);
      solver.solve (mp_equation, sol, rhs, m_preconditioner);
      deallog << "cg iters = " << solver_control.last_step () << std::endl;
    }
    //--------------------------------------------------------------------------
    // we use Richardson
    else if (solver_name == "richardson")
    {
      const unsigned int richardson_its = 1;
      this->richardson (sol, rhs, m_tol_i, richardson_its);
      //SolverRichardson<BlockVector<double> > solver (solver_control, m_vector_memory,
      //    SolverRichardson<BlockVector<double> >::AdditionalData (8));
      //solver.solve (mp_equation, sol, rhs, m_preconditioner);
      deallog << "richardson iters = " << richardson_its << std::endl;
    }
    //--------------------------------------------------------------------------
    /*SolverSelector<BlockVector<double> > solver();
     solver.select(solver_name);
     solver.set_control(solver_control);*/
    /*SolverSelector<BlockVector<double> >  solver(solver_name, solver_control, m_vector_memory);*/
    //--------------------------------------------------------------------------
    // We solve the system.
    //solver.solve (mp_equation, sol, rhs, m_preconditioner);
    //--------------------------------------------------------------------------
    //std::cout << "iterations = " << solver_control.last_step () << std::endl;
    m_total_wg_it += solver_control.last_step ();
  }

  template <int dim>
  void
  EqSolverWG<dim>::richardson (BlockVector<double> & sol,
      const BlockVector<double> & rhs,
      const double tol,
      const unsigned int max_it) const
  {
    BlockVector<double> m_sol_old (sol);
    for (unsigned int l = 0; l < max_it; ++l)
    {
      // Move sol_new to the container sol_old and start.
      m_sol_old.swap (sol);
      // equation.DL_invMS (sol, m_sol_old);
      mp_equation.vmult (sol, m_sol_old);
      sol.sadd (-1., m_sol_old);
      // Increase transport sweeps counter.
      ++m_total_wg_it;
      // New solution vector and its norm.
      sol += rhs;
      const double sol_norm = sol.l2_norm ();
      // Subtract the new solution to the old one to approximate the error.
      m_sol_old.sadd (-1, sol);
      const double err_phi = m_sol_old.l2_norm ();
      // If the error is small enough.. stop!.
      if (err_phi / sol_norm < tol)
      {
        break;
      }
    }
  }

  template class SolverEq<1> ;
  template class SolverEq<2> ;
  template class SolverEq<3> ;

  template class EqSolverWG<1> ;
  template class EqSolverWG<2> ;
  template class EqSolverWG<3> ;

} // end of namespace Forest
