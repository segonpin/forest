/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/neutronicmodule.cc
 * @brief  Implementation of class template NeutronicModule
 */

#include "neutronics/neutronicmodule.h"
#include "neutronics/hom_data.h"
#include "neutronics/homogenization.h"
#include "neutronics/state.h"

#include "algebra/eq_factory.h"
#include "algebra/equation.h"
#include "algebra/eigen_solver.h"

#include "input/input.h"
#include "input/input_settings.h"
#include "neutronics/extract_bv.h"
#include "output/forest_post_processing.h"
#include "algebra/eigen_prob_factory.h"

namespace Forest
{

  template <int dim>
  void
  NeutronicModule<dim>::run ()
  {
    // Preparing the (neutron transport) equation.
    const bool use_transport = mp_state.mp_data.mp_settings.use_transport ();
    std::shared_ptr<Equation<dim> > m_transport (
        EqFactory<dim>::New (mp_state, use_transport));

    // we use dsa to precondition transport. Otherwise does not make sense.
    bool use_dsa = mp_state.mp_data.mp_settings.use_dsa ();
    if (use_dsa)
    {
      deallog << "We are going to use DSA" << std::endl;
      // Preparing the (neutron diffusion) equation.
      std::shared_ptr<Equation<dim> > m_diffusion (
          EqFactory<dim>::New (mp_state, false));
      // Use diffusion as a preconditioner for transport?
      m_transport->set_prec (m_diffusion);
    }

    // Preparing the eigenvalue problem.
    std::shared_ptr<EigenProb<dim> > eigen_prob (
        EigenProbFactory<dim>::New (mp_state, m_transport));

    // Solving the neutron transport problem.
    EigenSolver<dim> eigen_solver (mp_state, eigen_prob);
    eigen_solver.solve_eigenproblem ();

    m_transport->memory_consumption();

    // Move this flag to the input
    bool homogenize = mp_state.mp_data.mp_settings.get_homogenize_flag ();
    if (homogenize and use_transport)
    {
      // Extract the boundary values for pins and assemblies.
      ExtractBV<dim> extract_bv (mp_state);
      // Homogenize and print values for generating discontinuity factors.
      Homogenization<dim> hom (mp_state);
      //  Postprocessing and printing
      PostProcessing<dim> output (mp_state);

      // Get the normalization constant and normalize interface data
      const double scaling = output.get_scaling ();
      extract_bv.scaling (scaling);
      extract_bv.scaling_moments (scaling);

      //  Printing the homogenized data
      HomData<dim> homdata (mp_state, hom, extract_bv);
    }
    else
    {
      //  Postprocessing and printing
      PostProcessing<dim> output (mp_state);
    }

    // Used to extract the boundary conditions and save them if required.
    /*if (mp_state.mp_data.mp_settings.get_print_ass_bcs_flag ())
     extract_bv.fill_bv (mp_state.m_scalar_flux);*/
  }

  template class NeutronicModule<1> ;
  template class NeutronicModule<2> ;
  template class NeutronicModule<3> ;

} // end of namespace Forest
