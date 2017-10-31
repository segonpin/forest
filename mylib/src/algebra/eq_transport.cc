/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eq_transport.cc
 * @brief  Implementation of class template EqTransport
 */

#include "algebra/eq_transport.h"

#include "utils/forest_utils_memory.h"      // for StatsMemory
#include "neutronics/trans_factory.h"       // for TransFactory
#include "neutronics/boundaryvalues.h"
#include "neutronics/boundary_values_factory.h"       // for TransFactory
#include "neutronics/fission_source_factory.h"  // for FissionSourceFactory
#include "neutronics/fission_spectrum_factory.h"  // for FissionSpectrumFactory
#include "neutronics/scattering_factory.h"  // for ScatteringFactory
#include "neutronics/state.h"               // for State
#include "neutronics/trans_op.h"            // for TransOp
#include "algebra/sn_vector.h"              // for SnVector
#include "angle/quadrature_base.h"
#include "angle/quadrature_factory.h"
#include "input/input.h"

#include "deal.II/base/thread_management.h" // for TaskGroup, new_task
#include "deal.II/base/exceptions.h"        // for Assert and ExcMessage

#include <vector>                           // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  EqTransport<dim>::EqTransport (State<dim> & state)
      : mp_state (state),
        m_solver (state, *this),
        m_matrix_free (mp_state.mp_data.mp_settings.get_matrix_free ()),
        m_ord (
            QuadratureFactory<dim>::build (
                mp_state.mp_data.mp_settings.get_sn_order (),
                mp_state.mp_data.mp_settings.get_quad_name ())),
        m_L0 (TransFactory<dim>::New (mp_state, m_ord, m_matrix_free)),
        m_BV (BoundaryValuesFactory<dim>::New (mp_state, m_ord, m_matrix_free)),
        m_Chi (FissionSpectrumFactory<dim>::New (mp_state, m_matrix_free)),
        m_XF (FissionSourceFactory<dim>::New (mp_state, m_matrix_free)),
        m_S (ScatteringFactory<dim>::New (mp_state, m_matrix_free)),
        m_M (m_ord, mp_state.get_n_leg_mom () - 1),
        m_D (m_ord, mp_state.get_n_leg_mom () - 1),
        m_g (0),
        m_external_preconditioner (false)
  {
    /* const unsigned int n_angles = mp_state.get_n_angles (); */
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();
    const unsigned int n_elem = mp_state.get_n_elem ();
    const unsigned int n_angles = m_ord->get_n_angles ();

    /* Temporary vectors. */
    m_v0.reinit (n_leg_mom, n_elem);
    m_v1.reinit (n_leg_mom, n_elem);
    /* Temporary vectors for single-thread. */
    m_v2.reinit (n_elem);
    m_v3.reinit (n_elem);
    /* Temporary vectors for multi-thread. */
    m_pv2.reinit (n_angles, n_elem);
    m_pv3.reinit (n_angles, n_elem);
    /* Temporary vectors for single-thread. */

    /* Memory used by the transport operator. */
    StatsMemory::instance ().add ("Transport matrix",
        m_L0->memory_consumption ()); /* Reporting memory usage. */
    /* Memory used by the fission operator. */
    StatsMemory::instance ().add ("Fission matrix",
        m_XF->memory_consumption ()); /* Reporting memory usage. */
    /* Memory used by the scattering operator. */
    StatsMemory::instance ().add ("Scattering matrix",
        m_S->memory_consumption ()); /* Reporting memory usage. */

    /* Set the size of the angular flux. */
    m_angular_sizes.clear ();
    m_angular_sizes.resize (n_angles, n_elem);

    /* Set the size of the leg moments and significant boundary dofs. */
    m_flux_bv_sizes.clear ();
    for (unsigned int m = 0; m < n_leg_mom; ++m)
    {
      m_flux_bv_sizes.push_back (n_elem);
    }
    for (unsigned int n = 0; n < n_angles; ++n)
    {
      m_flux_bv_sizes.push_back (m_BV->get_block_sizes (n));
    }

    /* Set the size of the fission density. */
    m_fission_density_sizes.clear ();
    m_fission_density_sizes.resize (1, n_elem);
  }

  template <int dim>
  bool
  EqTransport<dim>::no_upscatt (const unsigned int g) const
  {
    return m_S->no_upscatt (g);
  }

  template <int dim>
  const std::vector<unsigned int> &
  EqTransport<dim>::get_angular_sizes () const
  {
    return m_angular_sizes;
  }

  template <int dim>
  const std::vector<unsigned int> &
  EqTransport<dim>::get_flux_bv_sizes () const
  {
    return m_flux_bv_sizes;
  }

  template <int dim>
  const std::vector<unsigned int> &
  EqTransport<dim>::get_fission_density_sizes () const
  {
    return m_fission_density_sizes;
  }


  template <int dim>
  void
  EqTransport<dim>::vmult (SnVector & dst,
                           const SnVector & src) const
  {
    // are we working in parallel?
    //bool parallel = true;
    //if (parallel)
    //{
    //  this->ParDL_invMS (dst, src);
    //}
    //else
    //{
    this->DL_invMS (dst, src);
    //}

    /* Perform the U=V-U step. */
    /*
    for (unsigned int i = 0; i < dst.n_blocks (); ++i)
      if (dst.block (i).size () > 0)
        dst.block (i).sadd (-1., src.block (i));
    */
    dst.sadd (-1., 1.0, src);

  }



  template <int dim>
  void
  EqTransport<dim>::vmult (BlockVector<double> & dst,
                           const BlockVector<double> & src) const
  {
    /* are we working in parallel? */
    bool parallel = true;
    if (parallel)
    {
      this->ParDL_invMS (dst, src);
    }
    else
    {
      this->DL_invMS (dst, src);
    }

    /* Perform the U=V-U step. */
    /*
    for (unsigned int i = 0; i < dst.n_blocks (); ++i)
      if (dst.block (i).size () > 0)
        dst.block (i).sadd (-1., src.block (i));
    */
    dst.sadd (-1., src);

  }

  /** @brief Preconditioner ? */
  template <int dim>
  void
  EqTransport<dim>::prec (BlockVector<double> & dst,
          const BlockVector<double> & src) const
  {
    if (m_external_preconditioner)
    {
      std::cout << "EqTransport<dim>::prec; " << std::endl;
      /*
      std::cout << "src.n_blocks() = " << src.n_blocks () << std::endl;
      std::cout << "dst.n_blocks() = " << dst.n_blocks () << std::endl;
      */
      m_prec->set_group (m_g);
      m_S->in_scatt (m_v0, src, m_g);
      m_prec->solve (m_v1, m_v0, 1.e-5, 2000);
      dst = 0;
      dst.block (0) = m_v1.block (0);
      dst.sadd (+1.0, src);
    }
    else
    {
      //std::cout << "no preconditioner for transport" << std::endl;
      dst = src;
    }
  }

  /** @brief set the Preconditioner */
  template <int dim>
  void
  EqTransport<dim>::set_prec (std::shared_ptr<Equation<dim> > & prec) const
  {
    m_external_preconditioner = true;
    m_prec = prec;
  }

  template <int dim>
  void
  EqTransport<dim>::solve (BlockVector<double> & dst,
                           const BlockVector<double> & src,
                           const double tol_default,
                           const unsigned int max_it_default) const
  {
    /*std::cout << "EqTransport<dim>::solve; "<< std::endl;
     std::cout << "src.n_blocks() = " << src.n_blocks() << std::endl;
     std::cout << "dst.n_blocks() = " << dst.n_blocks() << std::endl;*/
    m_solver.solve (dst, src, tol_default, max_it_default);
  }

  /*template <int dim>
  void
  EqTransport<dim>::scatt_add (BlockVector<double> &dst,
                               const BlockVector<double> &src,
                               const unsigned int g,
                               const unsigned int h) const
  {
    // We directly apply the method from the scattering operator.
    m_S->vmult_add_g_h (dst, src, g, h);
  }*/

  template <int dim>
  void
  EqTransport<dim>::up_add (BlockVector<double> &dst,
                            const SnVector &src,
                            const unsigned int g) const
  {
    /* We directly apply the method from the scattering operator. */
    for (unsigned int h = g + 1; h < src.m_data.size (); ++h)
    {
      m_S->vmult_add_g_h (dst, src.m_data[h], g, h);
    }
  }

  template <int dim>
  void
  EqTransport<dim>::down_add (BlockVector<double> &dst,
                              const SnVector &src,
                              const unsigned int g) const
  {
    /* We directly apply the method from the scattering operator. */
    for (unsigned int h = 0; h < g; ++h)
    {
      m_S->vmult_add_g_h (dst, src.m_data[h], g, h);
    }
  }

  /** @brief access to the fission method. */
  /*template <int dim>
  void
  EqTransport<dim>::update_fission_source (SnVector &dst,
                                           const SnVector &src)
  {
    // m_XF.vmult (dst, src);
    const unsigned int n_elem = mp_state.get_n_elem ();
    SnVector tmp (1, 1, n_elem);
    this->calculate_fission_density (tmp, src);
    this->distribute_fission_rate (dst, tmp);
  }*/

  /** @brief access to the fission method. */
  template <int dim>
  void
  EqTransport<dim>::calculate_fission_density (SnVector &dst,
                                               const SnVector &src)
  {
    m_XF->vmult (dst, src);
  }

  /** @brief access to the fission method. */
  template <int dim>
  void
  EqTransport<dim>::distribute_fission_rate (SnVector &dst,
                                             const SnVector &src)
  {
    m_Chi->vmult (dst, src);
  }


  template <int dim>
  void
  EqTransport<dim>::DL_invMS (SnVector & dst,
      const SnVector &src) const
  {
    SnVector auxvec(src); // auxiliary vector with the right size
    //m_S->vmult(auxvec, src); // multiply and store the result in auxvec
    for (unsigned int g = 0; g < mp_state.get_n_groups (); g++)
    {
      m_S->in_scatt(auxvec, src, g);
      m_S->up_scatt_add(auxvec, src, g);
      m_S->down_scatt_add(auxvec, src, g);
      // i have to do this because the scattering overwrites the bcs
      unsigned int n_elem = mp_state.get_n_elem();
      unsigned int n_let_mom = mp_state.get_n_leg_mom();
      for (unsigned int i = n_elem*n_let_mom; i < src.m_data[g].size(); ++i)
        auxvec.m_data[g][i] = src.m_data[g][i];
    }
    DL_invM(dst, auxvec); // multiply by the transport sweep and store in dst
  }

  template <int dim>
  void
  EqTransport<dim>::DL_invMS (BlockVector<double> & dst,
                              const BlockVector<double> &src) const
  {
    // We need these constants.
    //const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    // Check if the energy group is out of range.
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));
    // if we have more than 1 legendre moments this is wrong
    Assert(mp_state.get_n_leg_mom () == 1, ExcMessage("n_leg_mom different from 1."));

    // multiply by in-scattering and put it in the temporary vector.
    BlockVector<double> tmp  = src; // I generate this vector to get the bcs
    m_S->in_scatt (tmp, src, m_g);

    // i have to do this because the scattering overwrites the bcs
    unsigned int n_elem = mp_state.get_n_elem();
    unsigned int n_let_mom = mp_state.get_n_leg_mom();
    for (unsigned int i = n_elem*n_let_mom; i < src.size(); ++i)
      tmp[i] = src[i];

    DL_invM (dst, tmp);
  }

  template <int dim>
  void
  EqTransport<dim>::DL_invM (SnVector & dst,
                             const SnVector &src) const
  {
    for (unsigned int g = 0; g < mp_state.get_n_groups(); g++)
    {
      m_g = g; // this is like set_group(g), but the last does not work
      DL_invM(dst.m_data[g], src.m_data[g]);
    }
  }


  template <int dim>
  void
  EqTransport<dim>::DL_invM (BlockVector<double> & dst,
                              const BlockVector<double> &src) const
  {
    /* We need these constants. */
    const unsigned int n_angles = m_ord->get_n_angles ();
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));
    /* if we have more than 1 legendre moments this is wrong */
    Assert(n_leg_mom == 1, ExcMessage("n_leg_mom different from 1."));

    /* Initialize destination vector to zero. */
    dst = 0;
    /* multiply by in-scattering and put it in the temporary vector. */
    m_v1 = src;
    for (unsigned int n = 0; n < n_angles; ++n)
    {
      /* Moments to direction and stored in m_v2. */
      m_M.vmult (m_v2, m_v1, n);
      /* Add reflective boundary conditions as incoming fluxes to m_v2. */
      m_BV->apply_reflective (m_v2, src, m_g, n, n_leg_mom);
      /* Transport sweep for direction n */
      m_L0->apply_L0_inv (m_v3, m_v2, m_g, n);
      /* Get outgoing fluxes for reflective boundary conditions. */
      m_BV->get_reflective (dst, m_v3, n, n_leg_mom);
      /* Apply directions-to-moments operator to recover the scalar flux. */
      m_D.vmult_add (dst, m_v3, n);
    }
  }

  template <int dim>
  void
  EqTransport<dim>::ParDL_invMS (BlockVector<double> & dst,
                                 const BlockVector<double> &src) const
  {
    /* We need these constants. */
    //const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));
    /* if we have more than 1 legendre moments this is wrong */
    Assert(mp_state.get_n_leg_mom () == 1, ExcMessage("n_leg_mom different from 1."));

    /* multiply by in-scattering and put it in the temporary vector. */
    //m_v1 = src; // I have to do this to copy the boundary conditions.
    BlockVector<double> tmp  = src; // I generate this vector to get the bcs

    /*
    unsigned int n_elem = mp_state.get_n_elem();
    unsigned int n_let_mom = mp_state.get_n_leg_mom();
    deallog << "tmp =";
    for (unsigned int i = n_elem*n_let_mom; i < tmp.size(); ++i)
      deallog << " " << tmp[i];
    deallog << std::endl;

    deallog << "src  =";
    for (unsigned int i = n_elem*n_let_mom; i < src.size(); ++i)
      deallog << " " << src[i];
    deallog << std::endl;
    */

    m_S->in_scatt (tmp, src, m_g);

    // i have to do this because the scattering overwrites the bcs
    unsigned int n_elem = mp_state.get_n_elem();
    unsigned int n_let_mom = mp_state.get_n_leg_mom();
    for (unsigned int i = n_elem*n_let_mom; i < src.size(); ++i)
      tmp[i] = src[i];

    ParDL_invM (dst, tmp) ;
  }

  template <int dim>
  void
  EqTransport<dim>::ParDL_invM (BlockVector<double> & dst,
                                const BlockVector<double> &src) const
  {
    /* We need these constants. */
    const unsigned int n_angles = m_ord->get_n_angles ();
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));
    /* if we have more than 1 legendre moments this is wrong */
    Assert(n_leg_mom == 1, ExcMessage("n_leg_mom different from 1."));

    /* Initialize destination vector to zero. */
    dst = 0;
    /* multiply by in-scattering and put it in the temporary vector. */

    m_v1 = src;

    const bool parallel = false;
    const bool parallel_sweep = true;

    /* Moments to direction and stored in m_v2. */
    if (parallel)
    {
      Threads::TaskGroup<void> task_group;
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        task_group += Threads::new_task (&MomentsToDirections<dim>::vmult, m_M,
            m_pv2.block (n), m_v1, n);
      }
      task_group.join_all ();
    }
    else
    {
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        m_M.vmult (m_pv2.block (n), m_v1, n);
      }
    }

    /* Add reflective boundary conditions as incoming fluxes to m_v2. */
    if (parallel)
    {
      BoundaryValues<dim> & BV = (*m_BV.get ());
      Threads::TaskGroup<void> task_group;
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        task_group += Threads::new_task (&BoundaryValues<dim>::apply_reflective,
            BV, m_pv2.block (n), src, m_g, n, n_leg_mom);
      }
      task_group.join_all ();
    }
    else
    {
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        m_BV->apply_reflective (m_pv2.block (n), src, m_g, n, n_leg_mom);
      }
    }

    /* Transport sweep for direction n */
    if (parallel_sweep)
    {
      TransOp<dim> & L0 = (*m_L0.get ());
      Threads::TaskGroup<void> task_group;
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        task_group += Threads::new_task (&TransOp<dim>::apply_L0_inv_g_iord, L0,
            m_pv3.block (n), m_pv2.block (n), m_g, n);
      }
      task_group.join_all ();
    }
    else
    {
      TransOp<dim> & L0 = (*m_L0.get ());
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        L0.apply_L0_inv (m_pv3.block (n), m_pv2.block (n), m_g, n);
      }
    }

    /* Get outgoing fluxes for reflective boundary conditions. */
    if (parallel)
    {
      BoundaryValues<dim> & BV = (*m_BV.get ());
      Threads::TaskGroup<void> task_group;
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        task_group += Threads::new_task (&BoundaryValues<dim>::get_reflective,
            BV, dst, m_pv3.block (n), n, n_leg_mom);
      }
      task_group.join_all ();
    }
    else
    {
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        m_BV->get_reflective (dst, m_pv3.block (n), n, n_leg_mom);
      }
    }

    /* Get outgoing fluxes for reflective boundary conditions. */
    if (parallel)
    {
      Threads::TaskGroup<void> task_group;
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        task_group += Threads::new_task (&DirectionsToMoments<dim>::vmult_add,
            m_D, dst, m_pv3.block (n), n);
      }
      task_group.join_all ();
    }
    else
    {
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        m_D.vmult_add (dst, m_pv3.block (n), n);
      }
    }

    /* Moments to direction and stored in m_v2. */
    /*for (unsigned int n = 0; n < n_angles; ++n)
     m_M.vmult (m_pv2.block(n), m_v1, n);*/

    /* Add reflective boundary conditions as incoming fluxes to m_v2. */
    /*for (unsigned int n = 0; n < n_angles; ++n)
     m_BV.apply_reflective (m_pv2.block(n), src, m_g, n, n_leg_mom);*/

    /* Transport sweep for direction n */
    /*for (unsigned int n = 0; n < n_angles; ++n)
     m_L0->apply_L0_inv (m_pv3.block(n), m_pv2.block(n), m_g, n);*/

    /* Get outgoing fluxes for reflective boundary conditions. */
    /*for (unsigned int n = 0; n < n_angles; ++n)
     m_BV.get_reflective (dst, m_pv3.block(n), n, n_leg_mom);*/

    /* Apply directions-to-moments operator to recover the scalar flux. */
    /*for (unsigned int n = 0; n < n_angles; ++n)
     m_D.vmult_add (dst, m_pv3.block(n), n);*/
  }

  /**
   * \f$ psi = L^{-1} M \left(S + \frac{1}{\lambda}*F \right)*\phi \f$
   */
  template <int dim>
  void
  EqTransport<dim>::get_angular (SnVector & dst,
                                 const SnVector & phi,
                                 const double lambda) const
  {
    /* We need these constants. */
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_angles = m_ord->get_n_angles ();
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();
    const unsigned int n_elem = mp_state.get_n_elem ();

    /* Initialize destination vector to zero. */
    dst = 0;

    /* Initialize source and set it equal to the scattering source. */
    SnVector source (n_groups, n_leg_mom, n_elem);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      m_S->in_scatt (source, phi, g);
      m_S->up_scatt_add (source, phi, g);
      m_S->down_scatt_add (source, phi, g);
    }

    /* Initialize fission density and set it equal to the fission source. */
    SnVector fission_density (n_groups, n_leg_mom, n_elem);
    //m_XF->vmult (fission_density, phi);

    SnVector tmp (1, 1, n_elem);
    m_XF->vmult (tmp, phi);
    m_Chi->vmult (fission_density, tmp);
    /* Add the fission source to the scattering source in source. */
    source.add (1 / lambda, fission_density);

    /* Perform sweeps for every group */
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      /* multiply by in-scattering and put it in the temporary vector. */
      for (unsigned int n = 0; n < n_angles; ++n)
      {
        /* Moments to direction and stored in m_v2. */
        m_M.vmult (m_v2, source.m_data[g], n);
        /* Add reflective boundary conditions as incoming fluxes to m_v2. */
        m_BV->apply_reflective (m_v2, phi.m_data[g], g, n, n_leg_mom);
        /* Transport sweep for direction n */
        m_L0->apply_L0_inv (dst.m_data[g].block (n), m_v2, g, n);
      }
    }
  }

  template <int dim>
  void
  EqTransport<dim>::set_rhs (BlockVector<double> & dst,
                             const BlockVector<double> &src) const
  {
    /* We need these constants. */
    const unsigned int n_angles = m_ord->get_n_angles ();
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();

    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));
    /* if we have more than 1 legendre moments this is wrong */
    Assert(n_leg_mom == 1, ExcMessage("n_leg_mom different from 1."));

    /* Initialize destination vector to zero. */
    dst = 0;
    /* First we prepare the right hand side */
    for (unsigned int n = 0; n < n_angles; ++n)
    {
      /* Moments to direction and stored in m_v2. */
      m_M.vmult (m_v2, src, n);
      /* Transport sweep for direction n */
      m_L0->apply_L0_inv (m_v3, m_v2, m_g, n);
      /* Get outgoing fluxes for reflective boundary conditions. */
      m_BV->get_reflective (dst, m_v3, n, n_leg_mom);
      /* Apply directions-to-moments operator to recover the scalar flux. */
      m_D.vmult_add (dst, m_v3, n);
    }
  }

  template <int dim>
  double
  EqTransport<dim>::memory_consumption () const
  {
    double cputime = m_L0->get_sweep_cputime();
    double walltime = m_L0->get_sweep_walltime();
    deallog << std::endl;
    deallog << "  CPU-time  for sweep is: " << cputime << std::endl;
    deallog << "  Wall-time for sweep is: " << walltime << std::endl << std::endl;

    double memory_consumption = 0;
    memory_consumption += 0;

    return memory_consumption;
  }

  template <int dim>
  unsigned int
  EqTransport<dim>::get_n_angles () const
  {
    return m_ord->get_n_angles ();
  }

  template <int dim>
  unsigned int
  EqTransport<dim>::get_n_groups() const
  {
    return mp_state.get_n_groups();
  }


  template class EqTransport<1> ;
  template class EqTransport<2> ;
  template class EqTransport<3> ;

} // end of namespace Forest
