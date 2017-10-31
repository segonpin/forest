/**
 * @brief  Example main program
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @date   Modified on 02/04/14
 */

#include "neutronics/manager.h"
#include "input/input.h"
#include "geometry/prob_geom.h"
#include "utils/forest_utils_logstream.h"

#include <neutronics/neutronicmodule.h> // for MultiPhysics
#include "neutronics/state.h"

#include <deal.II/base/exceptions.h>    // for Assert and ExcMessage
#include <deal.II/base/multithread_info.h>
#include "deal.II/base/mpi.h"

#include <iostream>
#include <string>

//-----------------------------------------------------------------------------
/**
 * @brief Template function for the setting and solution of the problem.
 * @param data
 */
template <int dim>
void
problem (Forest::Input & data)
{

  Forest::Manager<dim> manager (data);

  manager.set_lattice_problems();

  // Generating the geometry object
  Forest::ProbGeom<dim> geom (data);

  // Generating the state object
  Forest::State<dim> state (data, geom);

  // Generating the Neutronic Module
  Forest::NeutronicModule<dim> neutronics (state);
  neutronics.run ();
}

int
main (int argc,
      char **argv)
{
  //const unsigned int n_cores = dealii::MultithreadInfo::n_cores ();
  //std::cout << "n_cores = " << n_cores << std::endl;
  //If we do not limit the number or threads, it will be chosen by dealii 
  const unsigned int n_threads = dealii::MultithreadInfo::n_threads ();
  //std::cout << "n_threads = " << n_threads << std::endl;
  const unsigned int thread_limit = n_threads - 1;
  //dealii::MultithreadInfo::set_thread_limit (thread_limit);
  //dealii::MultithreadInfo::set_thread_limit (1);
  try
  {
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, thread_limit);
    const unsigned int depth_out = 2;
    // Preparing forest to be started.
    const std::string input_file = Forest::start_forest (argc, argv, depth_out);
    // All input data together.
    Forest::Input data (input_file);
    // assigning the spatial dimension of the problem.
    const unsigned int dim = data.mp_geom.get_dim ();
    // Choose the right dimension to solve the problem.
    switch (dim)
    {
      case 1:
        problem<1> (data);
        break;
      case 2:
        problem<2> (data);
        break;
      case 3:
        problem<3> (data);
        break;
      default:
        Assert(1 <= dim or dim <= 3, dealii::ExcMessage("Wrong dimension."))
        break;
    }
    // Report the timming and send the finishing message.
    Forest::finish_forest ();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl << exc.what ()
              << std::endl << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}
