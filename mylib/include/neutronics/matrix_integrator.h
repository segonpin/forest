/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/matrix_integrator.h
 * @brief  MatrixIntegrator class template declarations
 */

#ifndef MATRIX_INTEGRATOR_H
#define MATRIX_INTEGRATOR_H

#include "input/input.h"

#include "deal.II/meshworker/dof_info.h"
#include "deal.II/meshworker/integration_info.h"
#include "deal.II/meshworker/assembler.h"
#include "deal.II/meshworker/loop.h"

#include "deal.II/integrators/laplace.h"

#include "deal.II/base/config.h"
#include "deal.II/base/exceptions.h"
#include "deal.II/base/quadrature.h"
#include "deal.II/lac/full_matrix.h"
#include "deal.II/fe/mapping.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/meshworker/dof_info.h"

#include <memory>
#include <vector>

#ifndef DOXYGEN
namespace dealii { template <typename > class Vector; }
namespace dealii { template <int dim, int spacedim> class DoFHandler; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class MatrixIntegrator
   *
   * @ingroup ForestNeutronics
   *
   * @brief Matrix integrator for the discontinuous Galerkin
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
  class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:

    MatrixIntegrator (bool use_cell,
                      bool use_boundary,
                      bool use_face,
                      const unsigned int g,
                      const unsigned int h,
                      const Input &input,
                      bool only_mass = false)
        : MeshWorker::LocalIntegrator<dim, dim, double> (use_cell, use_boundary,
              use_face),
          m_g (g),
          m_h (h),
          mp_input (input),
          m_only_mass (only_mass)
    {
      /* Severification on the consistency*/
      Assert(g == h or only_mass,
          ExcMessage("Outside diagonal it is only mass matrix."));
    }

    unsigned int m_g;
    unsigned int m_h;
    const Input & mp_input;
    bool m_only_mass;

    void
    cell (MeshWorker::DoFInfo<dim> &dinfo,
          typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    boundary (MeshWorker::DoFInfo<dim> &dinfo,
              typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    face (MeshWorker::DoFInfo<dim> &dinfo1,
          MeshWorker::DoFInfo<dim> &dinfo2,
          typename MeshWorker::IntegrationInfo<dim> &info1,
          typename MeshWorker::IntegrationInfo<dim> &info2) const;

  private:

    void
    cell_matrix (FullMatrix<double> &M,
                 const FEValuesBase<dim> &fe,
                 const double diff) const;

    void
    nitsche_matrix (FullMatrix<double> &M,
                    const FEValuesBase<dim> &fe,
                    const double penalty,
                    const double diff) const;

    void
    robin_matrix (FullMatrix<double> &M,
                  const FEValuesBase<dim> &fe,
                  const double penalty,
                  const double diff,
                  const double factor) const;

    void
    ip_matrix (FullMatrix<double> &M11,
               FullMatrix<double> &M12,
               FullMatrix<double> &M21,
               FullMatrix<double> &M22,
               const FEValuesBase<dim> &fe1,
               const FEValuesBase<dim> &fe2,
               const double penalty,
               const double diff1,
               const double diff2) const;

    void
    mass_matrix (FullMatrix<double> &M,
                 const FEValuesBase<dim> &fe,
                 const double factor) const;

    double
    compute_penalty (const MeshWorker::DoFInfo<dim, dim, double> &dinfo1,
                     const MeshWorker::DoFInfo<dim, dim, double> &dinfo2,
                     const unsigned int deg1,
                     const unsigned int deg2,
                     const double coeff1,
                     const double coeff2) const;
  };

  template <int dim>
  void
  MatrixIntegrator<dim>::cell (MeshWorker::DoFInfo<dim> &dinfo,
                               typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    //std::cout << "cell only_mass = " << m_only_mass << std::endl;
    //std::cout << "m_g = " << m_g << "; m_h = " << m_h << std::endl<< std::endl;
    const unsigned int mat_id = dinfo.cell->material_id ();
    if (m_only_mass)
    {
      const double sigs = mp_input.mp_mat.m_xs.at (mat_id).m_sigmas[m_g][m_h];
      mass_matrix (dinfo.matrix (0, false).matrix, info.fe_values (), sigs);
    }
    else
    {
      const double diff = mp_input.mp_mat.m_xs.at (mat_id).m_diff[m_g];
      const double sigr = mp_input.mp_mat.m_xs.at (mat_id).m_sigmar[m_g];
      cell_matrix (dinfo.matrix (0, false).matrix, info.fe_values (), diff);
      mass_matrix (dinfo.matrix (0, false).matrix, info.fe_values (), sigr);
    }
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::boundary (MeshWorker::DoFInfo<dim> &dinfo,
                                   typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    //std::cout << "boundary " << std::endl;
    //std::cout << "m_g = " << m_g << "; m_h = " << m_h << std::endl<< std::endl;

    const unsigned int mat_id = dinfo.cell->material_id ();
    const double diff = mp_input.mp_mat.m_xs.at (mat_id).m_diff[m_g];
    const unsigned int deg = info.fe_values (0).get_fe ().tensor_degree ();

    const unsigned int fbi = dinfo.face->boundary_id ();
    const unsigned int b_ind = mp_input.mp_geom.m_core.m_bcs[fbi];
    // b_ind = 0 means zero incoming flux in transport, so albedo in diffusion
    // b_ind = 2 means reflective in transport, so zero current in diffusion
    if (b_ind != 2)
    {
      if (b_ind == 0)
      {
        const double albedo_val = 0.5;
        this->robin_matrix (dinfo.matrix (0, false).matrix, info.fe_values (),
            this->compute_penalty (dinfo, dinfo, deg, deg, diff, diff), diff,
            albedo_val);
      }
      else
      {
        this->nitsche_matrix (dinfo.matrix (0, false).matrix,
            info.fe_values (0),
            this->compute_penalty (dinfo, dinfo, deg, deg, diff, diff), diff);
      }
    }

  }

  template <int dim>
  void
  MatrixIntegrator<dim>::face (MeshWorker::DoFInfo<dim> &dinfo1,
                               MeshWorker::DoFInfo<dim> &dinfo2,
                               typename MeshWorker::IntegrationInfo<dim> &info1,
                               typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    //std::cout << "face " << std::endl;
    //std::cout << "m_g = " << m_g << "; m_h = " << m_h << std::endl<< std::endl;

    const unsigned int mat_id1 = dinfo1.cell->material_id ();
    const unsigned int mat_id2 = dinfo2.cell->material_id ();
    const double diff1 = mp_input.mp_mat.m_xs.at (mat_id1).m_diff[m_g];
    const double diff2 = mp_input.mp_mat.m_xs.at (mat_id2).m_diff[m_g];

    const unsigned int deg = info1.fe_values (0).get_fe ().tensor_degree ();

    this->ip_matrix (dinfo1.matrix (0, false).matrix,
        dinfo1.matrix (0, true).matrix, dinfo2.matrix (0, true).matrix,
        dinfo2.matrix (0, false).matrix, info1.fe_values (0),
        info2.fe_values (0),
        this->compute_penalty (dinfo1, dinfo2, deg, deg, diff1, diff2), diff1,
        diff2);
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::cell_matrix (FullMatrix<double> &M,
                                      const FEValuesBase<dim> &fe,
                                      const double diff) const
  {
    const unsigned int n_dofs = fe.dofs_per_cell;
    const unsigned int n_components = fe.get_fe ().n_components ();

    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      const double dx = fe.JxW (k) * diff;
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        for (unsigned int j = 0; j < n_dofs; ++j)
          for (unsigned int d = 0; d < n_components; ++d)
            M (i, j) += dx
                * (fe.shape_grad_component (j, k, d) * fe.shape_grad_component (
                       i, k, d));
      }
    }
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::nitsche_matrix (FullMatrix<double> &M,
                                         const FEValuesBase<dim> &fe,
                                         const double penalty,
                                         const double diff) const
  {
    const unsigned int n_dofs = fe.dofs_per_cell;
    const unsigned int n_comp = fe.get_fe ().n_components ();

    Assert(M.m() == n_dofs, ExcDimensionMismatch(M.m(), n_dofs));
    Assert(M.n() == n_dofs, ExcDimensionMismatch(M.n(), n_dofs));

    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      const double dx = fe.JxW (k);
      const Tensor<1, dim> &n = fe.normal_vector (k);
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          for (unsigned int d = 0; d < n_comp; ++d)
          {
            const double vi = fe.shape_value_component (i, k, d);
            const double dnvi = diff * n * fe.shape_grad_component (i, k, d);
            const double ui = fe.shape_value_component (j, k, d);
            const double dnui = diff * n * fe.shape_grad_component (j, k, d);
            M (i, j) += dx * (2. * penalty * vi * ui - dnvi * ui - dnui * vi);
          }
    }
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::robin_matrix (FullMatrix<double> &M,
                                       const FEValuesBase<dim> &fe,
                                       const double /*penalty*/,
                                       const double /*diff*/,
                                       const double factor) const
  {
    const unsigned int n_dofs = fe.dofs_per_cell;

    Assert(M.m() == n_dofs, ExcDimensionMismatch(M.m(), n_dofs));
    Assert(M.n() == n_dofs, ExcDimensionMismatch(M.n(), n_dofs));

    // this one has been working until now...
    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      const double dx = fe.JxW (k);
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const double vi = fe.shape_value (i, k);
          const double ui = fe.shape_value (j, k);
          M (i, j) += dx * factor * vi * ui;
        }
    }

    // lo que yo creo que debería de ser
    /*for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      const double dx = fe.JxW (k);
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const double vi = fe.shape_value (i, k);
          const double ui = fe.shape_value (j, k);
          for (unsigned int d = 0; d < fe.get_fe ().n_components (); ++d)
          {
            const Tensor<1, dim> &n = fe.normal_vector (k);
            const double dnvi = diff * n * fe.shape_grad_component (i, k, d);
            const double dnui = diff * n * fe.shape_grad_component (j, k, d);
            M (i, j) += dx * penalty
                        * (factor * vi * ui - dnvi * ui - dnui * vi)
                        + dx * (dnvi * ui + dnui * vi);
          }
        }
    }*/

  }

  template <int dim>
  void
  MatrixIntegrator<dim>::ip_matrix (FullMatrix<double> &M11,
                                    FullMatrix<double> &M12,
                                    FullMatrix<double> &M21,
                                    FullMatrix<double> &M22,
                                    const FEValuesBase<dim> &fe1,
                                    const FEValuesBase<dim> &fe2,
                                    const double penalty,
                                    const double diff1,
                                    const double diff2) const
  {
    const unsigned int n_dofs = fe1.dofs_per_cell;
    AssertDimension(M11.n(), n_dofs);
    AssertDimension(M11.m(), n_dofs);
    AssertDimension(M12.n(), n_dofs);
    AssertDimension(M12.m(), n_dofs);
    AssertDimension(M21.n(), n_dofs);
    AssertDimension(M21.m(), n_dofs);
    AssertDimension(M22.n(), n_dofs);
    AssertDimension(M22.m(), n_dofs);

    for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
    {
      const double dx = fe1.JxW (k);
      const Tensor<1, dim> &n = fe1.normal_vector (k);
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const double vi = fe1.shape_value (i, k);
          const double dnvi = diff1 * n * fe1.shape_grad (i, k);
          const double ve = fe2.shape_value (i, k);
          const double dnve = diff2 * n * fe2.shape_grad (i, k);
          const double ui = fe1.shape_value (j, k);
          const double dnui = diff1 * n * fe1.shape_grad (j, k);
          const double ue = fe2.shape_value (j, k);
          const double dnue = diff2 * n * fe2.shape_grad (j, k);

          M11 (i, j) += dx
              * (-.5 * dnvi * ui - .5 * dnui * vi + penalty * ui * vi);
          M12 (i, j) += dx
              * (.5 * dnvi * ue - .5 * dnue * vi - penalty * ue * vi);
          M21 (i, j) += dx
              * (-.5 * dnve * ui + .5 * dnui * ve - penalty * ui * ve);
          M22 (i, j) += dx
              * (.5 * dnve * ue + .5 * dnue * ve + penalty * ue * ve);
        }
      }
    }
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::mass_matrix (FullMatrix<double> &M,
                                      const FEValuesBase<dim> &fe,
                                      const double factor) const
  {
    const unsigned int n_dofs = fe.dofs_per_cell;
    const unsigned int n_components = fe.get_fe ().n_components ();
    for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
    {
      const double dx = fe.JxW (k) * factor;
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          for (unsigned int d = 0; d < n_components; ++d)
            M (i, j) += dx * fe.shape_value_component (j, k, d)
                        * fe.shape_value_component (i, k, d);
    }
  }

  template <int dim>
  double
  MatrixIntegrator<dim>::compute_penalty (const MeshWorker::DoFInfo<dim, dim,
                                              double> &dinfo1,
                                          const MeshWorker::DoFInfo<dim, dim,
                                              double> &dinfo2,
                                          const unsigned int deg1,
                                          const unsigned int deg2,
                                          const double coeff1,
                                          const double coeff2) const
  {
    const unsigned int normal1 =
        GeometryInfo<dim>::unit_normal_direction[dinfo1.face_number];
    const unsigned int normal2 =
        GeometryInfo<dim>::unit_normal_direction[dinfo2.face_number];
    const unsigned int deg1sq = (deg1 == 0) ? 1 : deg1 * (deg1 + 1);
    const unsigned int deg2sq = (deg2 == 0) ? 1 : deg2 * (deg2 + 1);

    double penalty1 = deg1sq / dinfo1.cell->extent_in_direction (normal1);
    double penalty2 = deg2sq / dinfo2.cell->extent_in_direction (normal2);
    if (dinfo1.cell->has_children () ^ dinfo2.cell->has_children ())
    {
      Assert(dinfo1.face == dinfo2.face, ExcInternalError());
      Assert(dinfo1.face->has_children(), ExcInternalError());
      penalty1 *= 2;
    }
    const double scal = 2.0;
    const double penalty = scal * 0.5 * (coeff1 * penalty1 + coeff2 * penalty2);
    return penalty;
  }

} // end of namespace Forest

#endif /* MATRIX_INTEGRATOR_H */
