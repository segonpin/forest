/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   algebra/eigen_prob_factory.cc
 * @brief  Implementation of class template EigenProbFactory
 */

#include "algebra/eigen_prob_factory.h"

#include "algebra/eigen_standard.h"
#include "algebra/eigen_fission_density.h"
#include "neutronics/state.h"
#include "input/input.h"

namespace Forest
{
  template <int dim>
  EigenProbPtr<dim>
  EigenProbFactory<dim>::New (State<dim> &state,
                              std::shared_ptr<Equation<dim> > &equation,
                              const bool default_use_fission_density)
  {
    bool use_fission_density = default_use_fission_density;
    const bool user_choice = state.mp_data.mp_settings.use_fission_density ();
    if (user_choice != use_fission_density)
    {

      use_fission_density = user_choice;
      std::cout << "User has chosen a different option from the suggested one"
                << std::endl;
    }

    if (use_fission_density)
    {
      return EigenProbPtr<dim> (new EigenFissionDensity<dim> (state, equation));
    }
    /* This is the default */
    return EigenProbPtr<dim> (new EigenStandard<dim> (state, equation));
  }

  template class EigenProbFactory<1> ;
  template class EigenProbFactory<2> ;
  template class EigenProbFactory<3> ;

} // end of namespace Forest

