/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/matrix_free_integrator.h
 * @brief  MatrixFreeIntegrator class template declarations
 */

#ifndef MATRIX_FREE_INTEGRATOR_H
#define MATRIX_FREE_INTEGRATOR_H

#include "neutronics/diff_op.h"
#include "neutronics/state_fwd.h"       // For class State&
#include "input/input.h"

#include "deal.II/lac/vector.h"          // for Vector
#include "deal.II/base/config.h"
#include "deal.II/base/exceptions.h"
#include "deal.II/base/quadrature.h"
#include "deal.II/fe/mapping.h"
#include "deal.II/fe/fe_values.h"

#include <memory>
#include <vector>

#ifndef DOXYGEN
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class MatrixFreeIntegrator
   *
   * @ingroup ForestNeutronics
   *
   * @brief Matrix Free integrator for the discontinuous Galerkin
   *
   * @details This class implements several functions:
   *
   *  - cell
   *  - boundary
   *  - face
   *
   * @todo update the documentation with weak formulation and matrix form
   */

  template <int dim>
  class MatrixFreeIntegrator
  {
  public:

    MatrixFreeIntegrator (const State<dim> & state,
                          const unsigned int g,
                          const unsigned int h)
        : mp_fe (state.get_fe ()),
          mp_mapping (state.get_mapping ()),
          mp_dof_handler (state.get_dof_handler ()),
          fe_degree (mp_fe.degree),
          quadrature (fe_degree + 1),
          update_flags (
              update_values | update_gradients | update_JxW_values
              | update_quadrature_points),
          fe_v (mp_mapping, mp_fe, quadrature, update_flags),
          face_quadrature (fe_degree + 1),
          face_update_flags (
              update_values | update_gradients | update_JxW_values
              | update_normal_vectors),
          fe_v_face (mp_mapping, mp_fe, face_quadrature, face_update_flags),
          fe_v_subface (mp_mapping, mp_fe, face_quadrature, face_update_flags),
          fe_v_subface_neighbor (mp_mapping, mp_fe, face_quadrature,
              face_update_flags),
          neighbor_update_flags (
              update_values | update_gradients | update_JxW_values
              | update_normal_vectors),
          fe_v_face_neighbor (mp_mapping, mp_fe, face_quadrature,
              neighbor_update_flags),
          dofs_per_cell (mp_fe.dofs_per_cell),
          m_g (g),
          m_h (h),
          mp_input (state.mp_data),
          m_only_mass (g != h)
    {
      cell_src.reinit (dofs_per_cell);
      cell_dst.reinit (dofs_per_cell);
      neighbor_src.reinit (dofs_per_cell);
      neighbor_dst.reinit (dofs_per_cell);
      local_dof_indices.resize (dofs_per_cell);
      neighbor_dof_indices.resize (dofs_per_cell);
    }

  private:

    // Short alias for the cell and face iterators
    typedef typename DoFHandler<dim, dim>::active_cell_iterator ac_iter;
    typedef typename DoFHandler<dim, dim>::cell_iterator c_iter;
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;

    const FiniteElement<dim, dim> & mp_fe;
    const Mapping<dim> & mp_mapping;
    const DoFHandler<dim, dim> & mp_dof_handler;
    unsigned int fe_degree;

    QGauss<dim> quadrature;
    UpdateFlags update_flags;
    mutable FEValues<dim> fe_v;

    QGauss<dim - 1> face_quadrature;
    UpdateFlags face_update_flags;
    mutable FEFaceValues<dim> fe_v_face;
    mutable FESubfaceValues<dim> fe_v_subface;
    mutable FESubfaceValues<dim> fe_v_subface_neighbor;

    UpdateFlags neighbor_update_flags;
    mutable FEFaceValues<dim> fe_v_face_neighbor;

    unsigned int dofs_per_cell;
    mutable Vector<double> cell_src, cell_dst, neighbor_src, neighbor_dst;
    mutable std::vector<unsigned int> local_dof_indices, neighbor_dof_indices;

    unsigned int m_g;
    unsigned int m_h;

    const Input & mp_input;
    bool m_only_mass;

    void
    cell (Vector<double> &dst,
          const Vector<double> &src,
          const FEValuesBase<dim> &fe,
          const c_iter &cell) const;

    void
    boundary (Vector<double> &dst,
              const Vector<double> &src,
              const FEValuesBase<dim> &fe,
              const c_iter &cell,
              const unsigned int face_no) const;
    void
    face (Vector<double> & cell_dst,
          Vector<double> & neighbor_dst,
          const Vector<double> & cell_src,
          const Vector<double> & neighbor_src,
          const FEValuesBase<dim> &fe1,
          const FEValuesBase<dim> &fe2,
          const c_iter &cell1,
          const unsigned int face1_no,
          const c_iter &cell2,
          const unsigned int face2_no) const;

  public:

    void
    vmult_add (Vector<double> &dst,
               const Vector<double> &src) const;

  private:

    void
    cell_matrix (Vector<double> &dst,
                 const Vector<double> &src,
                 const FEValuesBase<dim> &fe,
                 const double diff,
                 const double simgar) const;

    void
    nitsche_matrix (Vector<double> & dst,
                    const Vector<double> & src,
                    const FEValuesBase<dim> &fe,
                    const double penalty,
                    const double diff) const;

    void
    robin_matrix (Vector<double> & dst,
                  const Vector<double> & src,
                  const FEValuesBase<dim> &fe,
                  const double penalty,
                  const double diff,
                  const double factor) const;

    void
    ip_matrix (Vector<double> & cell_dst,
               Vector<double> & neighbor_dst,
               const Vector<double> & cell_src,
               const Vector<double> & neighbor_src,
               const FEValuesBase<dim> &fe1,
               const FEValuesBase<dim> &fe2,
               const double penalty,
               const double diff1,
               const double diff2) const;

    void
    mass_matrix (Vector<double> &dst,
                 const Vector<double> &src,
                 const FEValuesBase<dim> &fe,
                 const double factor) const;

    double
    compute_penalty (const c_iter &cell1,
                     const unsigned int face1_no,
                     const c_iter &cell2,
                     const unsigned int face2_no,
                     const unsigned int deg1,
                     const unsigned int deg2,
                     const double coeff1,
                     const double coeff2) const;
    double
    compute_penalty (const c_iter &cell,
                     const unsigned int face_no,
                     const unsigned int deg,
                     const double coeff) const;
  };

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::vmult_add (Vector<double> &dst,
                                        const Vector<double> &src) const
  {
    ac_iter cell = mp_dof_handler.begin_active (), endc = mp_dof_handler.end ();

    //std::cout << "MatrixFreeIntegrator" << std::endl;

    for (; cell != endc; ++cell)
    {
      // std::cout << "cell->user_index() = " << cell->user_index() << std::endl;
      // The local cell_dst equal to zero
      cell_dst = 0.;

      // Load the local cell_src
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        cell_src (i) = src[local_dof_indices[i]];

      // update the FeValues for this cell
      fe_v.reinit (cell);

      //----------------------------------------------------------------------
      // cell contribution
      //----------------------------------------------------------------------
      this->cell (cell_dst, cell_src, fe_v, cell);

      if (not m_only_mass)
      {
        // Then we run over all the neighbor cells
        for (unsigned int face_no = 0;
            face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          // Restart the FEFaceValues
          f_iter face = cell->face (face_no);
          fe_v_face.reinit (cell, face_no);

          //--------------------------------------------------------------------
          // boundary contribution
          //--------------------------------------------------------------------
          if (face->at_boundary ())
          {
            this->boundary (cell_dst, cell_src, fe_v_face, cell, face_no);
            continue;
          }
          //--------------------------------------------------------------------
          // we add each face share by two cells just one time.
          //--------------------------------------------------------------------
          c_iter neighbor = cell->neighbor (face_no);
          //--------------------------------------------------------------------
          // If face has children must be in the neighbor, not here
          //--------------------------------------------------------------------
          if (face->has_children ())
          {
            const unsigned int neighbor_face = cell->neighbor_face_no (face_no);
            for (unsigned int subface_no = 0;
                subface_no < face->number_of_children (); ++subface_no)
            {
              c_iter neigh_child = cell->neighbor_child_on_subface (face_no,
                  subface_no);
              if (neigh_child->user_index () < cell->user_index ())
              {
                continue;
              }
              fe_v_subface.reinit (cell, face_no, subface_no);
              fe_v_face_neighbor.reinit (neigh_child, neighbor_face);
              neigh_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                neighbor_src (i) = src[neighbor_dof_indices[i]];
              neighbor_dst = 0.;
              this->face (cell_dst, neighbor_dst, cell_src, neighbor_src,
                  fe_v_subface, fe_v_face_neighbor, cell, face_no, neigh_child,
                  neighbor_face);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                dst[neighbor_dof_indices[i]] += neighbor_dst (i);
            }
          }
          //--------------------------------------------------------------------
          // If face has children in this cell but not in the neighbor
          //--------------------------------------------------------------------
          else if (cell->neighbor_is_coarser (face_no))
          {
            if (neighbor->user_index () < cell->user_index ())
            {
              continue;
            }
            std::pair<unsigned int, unsigned int> neighbor_face_subface =
                cell->neighbor_of_coarser_neighbor (face_no);
            fe_v_face.reinit (cell, face_no);
            fe_v_subface.reinit (neighbor, neighbor_face_subface.first,
                neighbor_face_subface.second);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = src[neighbor_dof_indices[i]];
            neighbor_dst = 0.;
            this->face (cell_dst, neighbor_dst, cell_src, neighbor_src,
                fe_v_face, fe_v_subface, cell, face_no, neighbor,
                neighbor_face_subface.first);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              dst[neighbor_dof_indices[i]] += neighbor_dst (i);
          }
          //--------------------------------------------------------------------
          // If we are at the same level
          //--------------------------------------------------------------------
          else
          {
            if (neighbor->user_index () < cell->user_index ())
            {
              continue;
            }
            const unsigned int neighbor_face = cell->neighbor_of_neighbor (
                face_no);
            fe_v_face.reinit (cell, face_no);
            fe_v_face_neighbor.reinit (neighbor, neighbor_face);
            neighbor->get_dof_indices (neighbor_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              neighbor_src (i) = src[neighbor_dof_indices[i]];
            neighbor_dst = 0.;
            this->face (cell_dst, neighbor_dst, cell_src, neighbor_src,
                fe_v_face, fe_v_face_neighbor, cell, face_no, neighbor,
                neighbor_face);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              dst[neighbor_dof_indices[i]] += neighbor_dst (i);
          }
        }
      }
      // we map from local to global
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        dst[local_dof_indices[i]] += cell_dst (i);
    }
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::cell (Vector<double> &dst,
                                   const Vector<double> &src,
                                   const FEValuesBase<dim> &fe,
                                   const c_iter &cell) const
  {
    const unsigned int mat_id = cell->material_id ();
    if (m_only_mass)
    {
      const double sigs = mp_input.mp_mat.m_xs.at (mat_id).m_sigmas[m_g][m_h];
      mass_matrix (dst, src, fe, sigs);
    }
    else
    {
      const double diff = mp_input.mp_mat.m_xs.at (mat_id).m_diff[m_g];
      const double sigr = mp_input.mp_mat.m_xs.at (mat_id).m_sigmar[m_g];
      cell_matrix (dst, src, fe, diff, sigr);
    }
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::boundary (Vector<double> &dst,
                                       const Vector<double> &src,
                                       const FEValuesBase<dim> &fe,
                                       const c_iter &cell,
                                       const unsigned int face_no) const
  {
    const unsigned int mat_id = cell->material_id ();
    const double diff = mp_input.mp_mat.m_xs.at (mat_id).m_diff[m_g];
    const unsigned int deg = fe_degree;

    const unsigned int fbi = cell->face (face_no)->boundary_id ();
    const unsigned int b_ind = mp_input.mp_geom.m_core.m_bcs[fbi];
    // b_ind = 0 means zero incoming flux in transport, so albedo in diffusion
    // b_ind = 2 means reflective in transport, so zero current in diffusion
    if (b_ind != 2)
    {
      if (b_ind == 0)
      {
        const double albedo_val = 0.5;
        const double penalty = this->compute_penalty (cell, face_no, deg, diff);
        this->robin_matrix (dst, src, fe, penalty, diff, albedo_val);
      }
      else
      {
        const double penalty = this->compute_penalty (cell, face_no, deg, diff);
        this->nitsche_matrix (dst, src, fe, penalty, diff);
      }
    }

  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::face (Vector<double> & cell_dst,
                                   Vector<double> & neighbor_dst,
                                   const Vector<double> & cell_src,
                                   const Vector<double> & neighbor_src,
                                   const FEValuesBase<dim> &fe1,
                                   const FEValuesBase<dim> &fe2,
                                   const c_iter &cell1,
                                   const unsigned int face1_no,
                                   const c_iter &cell2,
                                   const unsigned int face2_no) const
  {
    const unsigned int mat_id1 = cell1->material_id ();
    const unsigned int mat_id2 = cell2->material_id ();
    const double diff1 = mp_input.mp_mat.m_xs.at (mat_id1).m_diff[m_g];
    const double diff2 = mp_input.mp_mat.m_xs.at (mat_id2).m_diff[m_g];
    const double penalty = this->compute_penalty (cell1, face1_no, cell2,
        face2_no, fe_degree, fe_degree, diff1, diff2);
    this->ip_matrix (cell_dst, neighbor_dst, cell_src, neighbor_src, fe1, fe2,
        penalty, diff1, diff2);
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::cell_matrix (Vector<double> &dst,
                                          const Vector<double> &src,
                                          const FEValuesBase<dim> &fe,
                                          const double diff,
                                          const double sigmar) const
  {
    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      // Prepare constants
      const double dx = fe.JxW (k);
      // Evaluate trial function
      Tensor<1, dim> dui;
      double ui = 0;
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        dui += fe.shape_grad (j, k) * src (j);
        ui += fe.shape_value (j, k) * src (j);
      }
      // Multiply by coefficient and quadrature weight
      dui *= dx * diff;
      ui *= dx * sigmar;
      // Multiply by test function
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        Tensor<1, dim> dvi = fe.shape_grad (i, k);
        double vi = fe.shape_value (i, k);
        dst (i) += dvi * dui + vi * ui;
      }
    }
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::nitsche_matrix (Vector<double> & dst,
                                             const Vector<double> & src,
                                             const FEValuesBase<dim> &fe,
                                             const double penalty,
                                             const double diff) const
  {
    for (unsigned int q = 0; q < fe_v_face.n_quadrature_points; ++q)
    {
      // Prepare constants
      const double dx = fe.JxW (q);
      const Tensor<1, dim> &n = fe.normal_vector (q);
      // Evaluate trial function
      double ui = 0;
      double dnui = 0;
      for (unsigned int j = 0; j < fe_v_face.dofs_per_cell; ++j)
      {
        ui += fe.shape_value (j, q) * cell_src (j);
        dnui += diff * n * fe.shape_grad (j, q) * src (j);
      }
      // Multiply by coefficient and quadrature weight
      ui *= dx;
      dnui *= dx;
      // Multiply by test function
      for (unsigned int i = 0; i < fe_v_face.dofs_per_cell; ++i)
      {
        const double vi = fe.shape_value (i, q);
        const double dnvi = diff * n * fe.shape_grad (i, q);
        dst (i) += 2. * penalty * vi * ui - dnvi * ui - dnui * vi;
      }
    }
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::robin_matrix (Vector<double> & dst,
                                           const Vector<double> & src,
                                           const FEValuesBase<dim> &fe,
                                           const double /* penalty */,
                                           const double /* diff */,
                                           const double factor) const
  {
    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      // Prepare constants
      const double dx = fe.JxW (k);
      // Evaluate trial function
      double ui = 0;
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        ui += fe.shape_value (j, k) * src (j);
      }
      // Multiply by coefficient and quadrature weight
      ui *= dx * factor;
      // Multiply by test function
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        double vi = fe.shape_value (i, k);
        dst (i) += vi * ui;
      }
    }

    // lo que yo creo que debería de ser
    /*for (unsigned int d = 0; d < fe.get_fe ().n_components (); ++d)
     {
     const Tensor<1, dim> &n = fe.normal_vector (k);
     const double dnvi = diff * n * fe.shape_grad_component (i, k, d);
     const double dnui = diff * n * fe.shape_grad_component (j, k, d);
     M (i, j) += dx * penalty
     * (factor * vi * ui - dnvi * ui - dnui * vi)
     + dx * (dnvi * ui + dnui * vi);
     }*/

  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::ip_matrix (Vector<double> & cell_dst,
                                        Vector<double> & neighbor_dst,
                                        const Vector<double> & cell_src,
                                        const Vector<double> & neighbor_src,
                                        const FEValuesBase<dim> &fe1,
                                        const FEValuesBase<dim> &fe2,
                                        const double penalty,
                                        const double diff1,
                                        const double diff2) const
  {
    for (unsigned int q = 0; q < fe1.n_quadrature_points; ++q)
    {
      // Prepare constants
      const double dx = fe1.JxW (q);
      const Tensor<1, dim> &n = fe1.normal_vector (q);
      // Evaluate trial function
      double ui = 0.;
      double dnui = 0.;
      double ue = 0.;
      double dnue = 0.;
      for (unsigned int j = 0; j < fe1.dofs_per_cell; ++j)
      {
        ui += fe1.shape_value (j, q) * cell_src (j);
        dnui += n * fe1.shape_grad (j, q) * cell_src (j);
        ue += fe2.shape_value (j, q) * neighbor_src (j);
        dnue += n * fe2.shape_grad (j, q) * neighbor_src (j);
      }
      // Multiply by coefficient and quadrature weight
      ui *= dx;
      dnui *= diff1 * dx;
      ue *= dx;
      dnue *= diff2 * dx;
      // Multiply by test function
      for (unsigned int i = 0; i < fe1.dofs_per_cell; ++i)
      {
        const double vi = fe1.shape_value (i, q);
        const double dnvi = diff1 * n * fe1.shape_grad (i, q);
        const double ve = fe2.shape_value (i, q);
        const double dnve = diff2 * n * fe2.shape_grad (i, q);
        cell_dst (i) += -.5 * dnvi * ui - .5 * dnui * vi + penalty * ui * vi;
        cell_dst (i) += .5 * dnvi * ue - .5 * vi * dnue - penalty * vi * ue;
        neighbor_dst (i) += -.5 * dnve * ui + .5 * ve * dnui
            - penalty * ve * ui;
        neighbor_dst (i) += .5 * dnve * ue + .5 * ve * dnue + penalty * ve * ue;
      }
    }
  }

  template <int dim>
  void
  MatrixFreeIntegrator<dim>::mass_matrix (Vector<double> &dst,
                                          const Vector<double> &src,
                                          const FEValuesBase<dim> &fe,
                                          const double factor) const
  {
    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      // Prepare constants
      const double dx = fe.JxW (k);
      // Evaluate trial function
      double ui = 0;
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        ui += fe.shape_value (j, k) * src (j);
      }
      // Multiply by coefficient and quadrature weight
      ui *= dx * factor;
      // Multiply by test function
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        double vi = fe.shape_value (i, k);
        dst (i) += vi * ui;
      }
    }
  }

  template <int dim>
  double
  MatrixFreeIntegrator<dim>::compute_penalty (const c_iter &cell1,
                                              const unsigned int face1_no,
                                              const c_iter &cell2,
                                              const unsigned int face2_no,
                                              const unsigned int deg1,
                                              const unsigned int deg2,
                                              const double coeff1,
                                              const double coeff2) const
  {
    const unsigned int normal1 =
        GeometryInfo<dim>::unit_normal_direction[face1_no];
    const unsigned int normal2 =
        GeometryInfo<dim>::unit_normal_direction[face2_no];
    const unsigned int deg1sq = (deg1 == 0) ? 1 : deg1 * (deg1 + 1);
    const unsigned int deg2sq = (deg2 == 0) ? 1 : deg2 * (deg2 + 1);

    double penalty1 = deg1sq / cell1->extent_in_direction (normal1);
    double penalty2 = deg2sq / cell2->extent_in_direction (normal2);
    if (cell1->has_children () ^ cell2->has_children ())
    {
      Assert(cell1->face(face1_no) == cell2->face(face2_no), ExcInternalError());
      Assert(cell1->face(face1_no)->has_children(), ExcInternalError());
      penalty1 *= 2;
    }
    const double scal = 2.0;
    const double penalty = scal * 0.5 * (coeff1 * penalty1 + coeff2 * penalty2);
    return penalty;
  }

  template <int dim>
  double
  MatrixFreeIntegrator<dim>::compute_penalty (const c_iter &cell,
                                              const unsigned int face_no,
                                              const unsigned int deg,
                                              const double coeff) const
  {
    const unsigned int normal =
        GeometryInfo<dim>::unit_normal_direction[face_no];
    const unsigned int degsq = (deg == 0) ? 1 : deg * (deg + 1);
    const double scal = 2.0;
    double penalty = scal * degsq * coeff / cell->extent_in_direction (normal);
    return penalty;
  }

} // end of namespace Forest

#endif /* MATRIX_FREE_INTEGRATOR_H */
