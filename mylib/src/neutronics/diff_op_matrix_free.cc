/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/diff_op_matrix_free.cc
 * @brief  Implementation of class template DiffOpMatrixFree
 */

#include "neutronics/diff_op_matrix_free.h"

#include "neutronics/matrix_free_integrator.h"
#include "utils/forest_utils_dealii.h"
#include "geometry/path_to_gridgenerator_dealii.h"
#include "geometry/prob_geom.h"
#include "input/input.h"
#include "neutronics/state.h"

#include "deal.II/base/config.h"          // for DEAL_II_VERSION_GTE
#include "deal.II/base/geometry_info.h"  // for GeometryInfo
#include "deal.II/base/point.h"          // for Point
#include "deal.II/base/tensor.h"         // for Tensor
#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/dofs/dof_tools.h"      // for FiniteElement, Mapping
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc
#include "deal.II/lac/vector.h"          // for Vector
#include "deal.II/lac/full_matrix.h"     // for FullMatrix
#include "deal.II/fe/fe_values.h"        // for FEFaceValues, etc

#include "deal.II/meshworker/dof_info.h"
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/assembler.h"
#include "deal.II/meshworker/loop.h"
#include "deal.II/integrators/laplace.h"
#include "deal.II/lac/sparse_matrix.h"

#include <utility>                       // for pair
#include <vector>                        // for vector

namespace Forest
{
  using namespace dealii;

  template <int dim>
  DiffOpMatrixFree<dim>::DiffOpMatrixFree (State<dim> & state)
      : mp_state (state),
        mp_fe (state.get_fe ()),
        mp_mapping (state.get_mapping ()),
        mp_dof_handler (mp_state.get_dof_handler ())
  {
    this->setup_system ();
    this->assemble_system ();
    this->setup_prec ();
    this->assemble_prec ();
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::vmult (BlockVector<double> & phi,
                                const BlockVector<double> & source,
                                const unsigned int g)
  {
    //matfree_integrator (
    //    new MatrixFreeIntegrator<dim> (mp_state, g, g));
    matfree_integrator[g]->vmult_add (phi.block (0), source.block (0));
    //this->vmult_matrix_free_g_h (phi, source, g, h);
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::prec (BlockVector<double> & phi,
                               const BlockVector<double> & source,
                               const unsigned int /* g */) const
  {
    phi = source;
    //const unsigned int n_m = 1;
    //for (unsigned int i = 0; i < mp_dof_handler.n_dofs (); ++i)
    //  phi.block (0)[i] = source.block (0)[i] * m_prec_matrix[g].block (0)[i];
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::vmult_matrix_free_g_h (BlockVector<double> & dst,
                                                const BlockVector<double> & src,
                                                const unsigned int g)
  {
    matfree_integrator[g]->vmult_add (dst.block (0), src.block (0));
    //deallog << "m_Diff_matrix[g][h].vmult_add (dst, src);" << std::endl;
    //m_Diff_matrix[g].vmult_add (dst, src);

  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::setup_system ()
  {
    const unsigned int n_g = mp_state.get_n_groups ();
    for (unsigned int g = 0; g < n_g; ++g)
    {
      matfree_integrator.push_back (
        std::shared_ptr<MatrixFreeIntegrator<dim> > (
          new MatrixFreeIntegrator<dim> (mp_state, g, g)));
    }
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::assemble_system ()
  {
    // do nothing
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::setup_prec ()
  {
    // do nothing
  }

  template <int dim>
  void
  DiffOpMatrixFree<dim>::assemble_prec ()
  {
    // do nothing
  }

  template <int dim>
  unsigned int
  DiffOpMatrixFree<dim>::get_n_groups () const
  {
    return mp_state.get_n_groups ();
  }

  template <int dim>
  double
  DiffOpMatrixFree<dim>::memory_consumption () const
  {
    //const unsigned int n_groups = mp_state.get_n_groups ();
    double memory_consumption = 0;
    return memory_consumption;
  }

  template class DiffOpMatrixFree<1> ;
  template class DiffOpMatrixFree<2> ;
  template class DiffOpMatrixFree<3> ;

} // end of namespace Forest
