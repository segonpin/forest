/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/trans_bc.cc
 * @brief  Implementation of class template TransBc
 */

#include "neutronics/trans_bc.h"

#include "angle/quadrature_base.h"
#include "neutronics/state.h"
#include "input/input.h"
#include "input/input_geom.h"
#include "algebra/sn_vector.h"             // for SnVector
#include "utils/forest_utils_dealii.h"  // for vector_to_point

#include "deal.II/base/config.h"        // for DEAL_II_VERSION_GTE
#include "deal.II/base/table.h"         // for Table
#include "deal.II/fe/fe_values.h"       // for FEFaceValues, FEValues
#include "deal.II/lac/full_matrix.h"    // for FullMatrix
#include "deal.II/base/geometry_info.h" // for GeometryInfo
#include "deal.II/base/point.h"         // for Point
#if DEAL_II_VERSION_GTE(8,4,0)
#include "deal.II/base/tensor.h"        // for Tensor
#endif
#include "deal.II/base/quadrature_lib.h"// for QGauss
#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/dofs/dof_tools.h"     // for FiniteElement
#include "deal.II/fe/fe_update_flags.h" // for operator|, UpdateFlags, etc
#include "deal.II/lac/block_vector.h"   // for BlockVector
#include "deal.II/lac/vector.h"         // for Vector

#include <vector>
#include <memory>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  TransBc<dim>::TransBc (State<dim> & state,
      std::shared_ptr<QuadratureBase<dim> > & ord)
      : mp_state (state),
        mp_ord (ord),
        mp_dof_handler (state.get_dof_handler ())
  {
    m_n_angles = mp_ord->get_n_angles();
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR (SnVector & source, const SnVector & phi_bc)
  {
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      apply_LR_g (source.m_data[g], phi_bc.m_data[g], g);
    }
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR (BlockVector<double> & source, const SnVector & phi_bc,
                          const unsigned int g)
  {
    apply_LR_g (source, phi_bc.m_data[g], g);
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR (Vector<double> & source, const SnVector & phi_bc,
                          const unsigned int g, const unsigned int i_ord)
  {
    apply_LR_g_iord (source, phi_bc.m_data[g], g, i_ord);
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR (Vector<double> & source,
                          const BlockVector<double> & phi_bc,
                          const unsigned int g, const unsigned int i_ord)
  {
    apply_LR_g_iord (source, phi_bc, g, i_ord);
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR_g (BlockVector<double> & source,
                            const BlockVector<double> & sol_old,
                            const unsigned int g)
  {
    // we should enter this loop only if we have specular reflection
    // boundary conditions.
    for (unsigned int i_ord = 0; i_ord < m_n_angles; ++i_ord)
    {
      apply_LR_g_iord (source.block (i_ord), sol_old, g, i_ord);
    }
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR_g_iord (Vector<double> & source,
                                 const BlockVector<double> & sol_old,
                                 const unsigned int g, const unsigned int i_ord)
  {
    // we should enter this loop only if we have specular reflection
    // boundary conditions.
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
    {
      unsigned int j_ord = mp_ord->get_out_ind (f, i_ord);
      if (mp_state.mp_data.mp_geom.m_core.m_bcs[f] == 2 and
          j_ord != (unsigned int) -1)
      {
        apply_LR_g_iord_f (source, sol_old.block (j_ord), g, i_ord, f);
      }
    }
  }

  template <int dim>
  void
  TransBc<dim>::apply_LR_g_iord_f (Vector<double> & dst,
                                   const Vector<double> & src,
                                   const unsigned int /* g */,
                                   const unsigned int i_ord,
                                   const unsigned int boundary_face)
  {
    const Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));

    // we get the finite elements
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());

    QGauss<dim> quadrature (fe.degree + 1);
    const UpdateFlags update_flags = update_values | update_JxW_values
                                     | update_quadrature_points;
    FEValues<dim> fe_values (mp_state.get_mapping (), fe, quadrature,
        update_flags);

    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
                                          | update_quadrature_points
                                          | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mp_state.get_mapping (), fe,
        face_quadrature, face_update_flags);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_src (dofs_per_cell), cell_dst (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim, dim>::active_cell_iterator cell =
        mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();
    for (; cell != endc; ++cell)
    {
      local_matrix = 0;

      bool apply_local_matrix = false;
      for (unsigned int face_no = 0;
          face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        typename DoFHandler<dim, dim>::face_iterator face = cell->face (
            face_no);
        unsigned int fbi = face->boundary_id ();
        if (face->at_boundary () and fbi == boundary_face)
        {
          apply_local_matrix = true;
          fe_values_face.reinit (cell, face_no);
          const std::vector<double> &JxW = fe_values_face.get_JxW_values ();
#if DEAL_II_VERSION_GTE(8,4,0)
          const std::vector<Tensor<1,dim> > &normals = 
            fe_values_face.get_all_normal_vectors ();
#else
          const std::vector<Point<dim> > &normals =
              fe_values_face.get_normal_vectors ();
#endif

          for (unsigned int q_point = 0;
              q_point < fe_values_face.n_quadrature_points; ++q_point)
          {
            const double beta_n = beta * normals[q_point];
            if (beta_n < 0)
            {
              for (unsigned int i = 0; i < fe_values_face.dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < fe_values_face.dofs_per_cell; ++j)
                {
                  local_matrix (i, j) -= beta_n
                      * fe_values_face.shape_value (j, q_point)
                      * fe_values_face.shape_value (i, q_point) * JxW[q_point];
                }
              }
            }
          }
        }
      }
      // we only perform the multiplication if necessary
      if (apply_local_matrix)
      {
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          cell_src (i) = src[local_dof_indices[i]];
        }
        local_matrix.vmult (cell_dst, cell_src);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          dst[local_dof_indices[i]] += cell_dst (i);
        }
      }
    }

  }

  template class TransBc<1> ;
  template class TransBc<2> ;
  template class TransBc<3> ;

} // end of namespace Forest
