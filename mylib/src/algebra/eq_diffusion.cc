/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eq_diffusion.cc
 * @brief  Implementation of class template EqDiffusion
 */

#include "algebra/eq_diffusion.h"

#include "utils/forest_utils_memory.h"      // for StatsMemory
#include "neutronics/diff_factory.h"       // for DiffFactory
#include "neutronics/diff_op.h"            // for DiffOp
#include "neutronics/fission_source_factory.h"  // for FissionSourceFactory
#include "neutronics/fission_spectrum_factory.h"  // for FissionSpectrumFactory
#include "neutronics/scattering_factory.h"  // for ScatteringFactory
#include "neutronics/state.h"               // for State
#include "algebra/sn_vector.h"              // for SnVector
#include "input/input.h"

#include "deal.II/base/thread_management.h" // for TaskGroup, new_task
#include "deal.II/base/exceptions.h"        // for Assert and ExcMessage

#include <vector>                           // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  EqDiffusion<dim>::EqDiffusion (State<dim> & state)
      : mp_state (state),
        m_solver (state, *this),
        //        m_matrix_free (mp_state.mp_data.mp_settings.get_matrix_free ()),
        m_matrix_free (false),
        //m_Diff (DiffFactory<dim>::New (mp_state, m_matrix_free)),
        m_Diff (DiffFactory<dim>::New (mp_state, false)),
        m_Chi (FissionSpectrumFactory<dim>::New (mp_state, m_matrix_free)),
        m_XF (FissionSourceFactory<dim>::New (mp_state, m_matrix_free)),
        m_S (ScatteringFactory<dim>::New (mp_state, m_matrix_free)),
        m_g (0),
        m_external_preconditioner (false)
  {

    /* const unsigned int n_angles = mp_state.get_n_angles (); */
    const unsigned int n_leg_mom = mp_state.get_n_leg_mom ();
    const unsigned int n_elem = mp_state.get_n_elem ();

    /* Temporary vectors. */
    m_v1.reinit (n_leg_mom, n_elem);
    /* Temporary vectors for single-thread. */
    m_v2.reinit (n_elem);
    m_v3.reinit (n_elem);

    /* Memory used by the transport operator. */
    StatsMemory::instance ().add ("Diffusion matrix",
        m_Diff->memory_consumption ()); /* Reporting memory usage. */
    /* Memory used by the fission operator. */
    StatsMemory::instance ().add ("Fission matrix",
        m_XF->memory_consumption ()); /* Reporting memory usage. */
    /* Memory used by the scattering operator. */
    StatsMemory::instance ().add ("Scattering matrix",
        m_S->memory_consumption ()); /* Reporting memory usage. */

    /* Set the size of the angular flux. */
    m_angular_sizes.clear ();
    m_angular_sizes.resize (1, n_elem);

    /* Set the size of the leg moments and significant boundary dofs. */
    m_flux_bv_sizes.clear ();
    for (unsigned int m = 0; m < n_leg_mom; ++m)
    {
      m_flux_bv_sizes.push_back (n_elem);
    }
    // for (unsigned int n = 0; n < n_angles; ++n)
    //   m_flux_bv_sizes.push_back (m_BV->get_block_sizes (n));

    /* Set the size of the fission density. */
    m_fission_density_sizes.clear ();
    m_fission_density_sizes.resize (1, n_elem);
  }

  template <int dim>
  bool
  EqDiffusion<dim>::no_upscatt (const unsigned int g) const
  {
    return m_S->no_upscatt (g);
  }

  template <int dim>
  const std::vector<unsigned int> &
  EqDiffusion<dim>::get_angular_sizes () const
  {
    return m_angular_sizes;
  }

  template <int dim>
  const std::vector<unsigned int> &
  EqDiffusion<dim>::get_flux_bv_sizes () const
  {
    return m_flux_bv_sizes;
  }
  template <int dim>
  const std::vector<unsigned int> &
  EqDiffusion<dim>::get_fission_density_sizes () const
  {
    return m_fission_density_sizes;
  }

  template <int dim>
  void
  EqDiffusion<dim>::vmult (BlockVector<double> & dst,
                           const BlockVector<double> & src) const
  {
    /*std::cout << "EqDiffusion<dim>::vmult; "<< std::endl;
    std::cout << "src.n_blocks() = " << src.n_blocks() << std::endl;
    std::cout << "dst.n_blocks() = " << dst.n_blocks() << std::endl;*/

    /* are we working in parallel? */
    bool parallel = true;
    if (parallel)
    {
      this->ParDiffS (dst, src);
    }
    else
    {
      this->DiffS (dst, src);
    }
  }

  /** @brief Preconditioner ? */
  template <int dim>
  void
  EqDiffusion<dim>::prec (BlockVector<double> & dst,
                          const BlockVector<double> & src) const
  {
    if (m_external_preconditioner)
    {
      //std::cout << "external preconditioner for diffusion" << std::endl;
      m_prec->set_group (m_g);
      m_prec->solve (dst, src);
    }
    else
    {
      //std::cout << "default preconditioner for diffusion" << std::endl;
      /*std::cout << "EqDiffusion<dim>::prec; "<< std::endl;
      std::cout << "src.n_blocks() = " << src.n_blocks() << std::endl;
      std::cout << "dst.n_blocks() = " << dst.n_blocks() << std::endl;*/
      m_Diff->prec (dst, src, m_g);
    }
  }

  /** @brief set the Preconditioner */
  template <int dim>
  void
  EqDiffusion<dim>::set_prec (std::shared_ptr<Equation<dim> > & prec) const
  {
    m_external_preconditioner = true;
    m_prec = prec;
  }

  template <int dim>
  void
  EqDiffusion<dim>::solve (BlockVector<double> & dst,
                           const BlockVector<double> & src,
                           const double tol_default,
                           const unsigned int max_it_default) const
  {
    /*std::cout << "EqDiffusion<dim>::solve; "<< std::endl;
    std::cout << "src.n_blocks() = " << src.n_blocks() << std::endl;
    std::cout << "dst.n_blocks() = " << dst.n_blocks() << std::endl;*/
    m_solver.solve (dst, src, tol_default, max_it_default);
  }

  /** @todo do it with vmult_add? */
  template <int dim>
  void
  EqDiffusion<dim>::up_add (BlockVector<double> &dst,
                            const SnVector &src,
                            const unsigned int g) const
  {
    /* We directly apply the method from the scattering operator. */
    for (unsigned int h = g + 1; h < src.m_data.size (); ++h)
    {
      m_S->vmult_add_g_h (dst, src.m_data[h], g, h);
    }
  }

  /** @todo do it with vmult_add? */
  template <int dim>
  void
  EqDiffusion<dim>::down_add (BlockVector<double> &dst,
                              const SnVector &src,
                              const unsigned int g) const
  {
    /* We directly apply the method from the scattering operator. */
    for (unsigned int h = 0; h < g; ++h)
    {
      m_S->vmult_add_g_h (dst, src.m_data[h], g, h);
    }
  }

  /** @todo Remove the auxiliary vector from here.. */
  /*template <int dim>
  void
  EqDiffusion<dim>::update_fission_source (SnVector &dst,
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
  EqDiffusion<dim>::calculate_fission_density (SnVector &dst,
                                            const SnVector &src)
  {
    m_XF->vmult (dst, src);
  }

  /** @brief access to the fission method. */
  template <int dim>
  void
  EqDiffusion<dim>::distribute_fission_rate (SnVector &dst,
                                             const SnVector &src)
  {
    m_Chi->vmult (dst, src);
  }

  template <int dim>
  void
  EqDiffusion<dim>::DiffS (BlockVector<double> & dst,
                           const BlockVector<double> &src) const
  {
    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));

    /* Initialize destination vector to zero. */
    dst = 0;
    /* Transport sweep for direction n */
    m_Diff->vmult (dst, src, m_g);
  }

  /** @todo Parallelization in angle makes no sense here.
   * Maybe other parallelization*/
  template <int dim>
  void
  EqDiffusion<dim>::ParDiffS (BlockVector<double> & dst,
                              const BlockVector<double> &src) const
  {
    this->DiffS (dst, src);
  }

  /** @todo can we add some kind of approximation? */
  template <int dim>
  void
  EqDiffusion<dim>::get_angular (SnVector & dst,
                                 const SnVector & phi,
                                 const double /* lambda */) const
  {
    //Assert(false, ExcMessage("Not implemented."));
    dst = phi;
  }

  template <int dim>
  void
  EqDiffusion<dim>::set_rhs (BlockVector<double> & dst,
                             const BlockVector<double> & src) const
  {
    /* Check if the energy group is out of range. */
    Assert(m_g < mp_state.get_n_groups (),
        ExcMessage("group index out of range."));

    /* We just copy the vector */
    dst = src;
  }

  template <int dim>
  double
  EqDiffusion<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption += 0;
    return memory_consumption;
  }

  template <int dim>
  unsigned int
  EqDiffusion<dim>::get_n_angles () const
  {
    //Assert(false, ExcMessage("Diffusion equation has no angles."));
    return 1;
  }

  template <int dim>
  unsigned int
  EqDiffusion<dim>::get_n_groups() const
  {
    return mp_state.get_n_groups();
  }


  template class EqDiffusion<1> ;
  template class EqDiffusion<2> ;
  template class EqDiffusion<3> ;

} // end of namespace Forest
