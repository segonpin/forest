/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/boundaryvalues.cc
 * @brief  Implementation of class template BoundaryValues
 */

#include "neutronics/boundaryvalues.h"

#include "angle/quadrature_base.h"
#include "utils/forest_utils_dealii.h"
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "input/input.h"
#include "neutronics/state.h"
#include "geometry/prob_geom.h"

#include "deal.II/base/tensor.h"        // for Tensor
#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/fe/fe_values.h"        // for FEValues
#include "deal.II/fe/fe.h"               // for FiniteElement
#include "deal.II/base/geometry_info.h"  // for GeometryInfo
#include "deal.II/base/point.h"          // for Point
#include "deal.II/base/types.h"          // for global_dof_index
#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc
#include "deal.II/fe/mapping_q1.h"       // for MappingQ1
#include "deal.II/lac/vector.h"           // for Vector
#include "deal.II/lac/block_vector.h"    // for BlockVector
#include "deal.II/base/exceptions.h"     // for Assert and ExcMessage

#include <utility>                       // for pair
#include <map>

namespace Forest
{

  using namespace dealii;

  template <int dim>
  BoundaryValues<dim>::BoundaryValues (State<dim> & state,
                                       std::shared_ptr<QuadratureBase<dim> > & ord)
      : mp_state (state),
        mp_ord (ord),
        m_fe_sup_on_face (GeometryInfo<dim>::faces_per_cell)
  {
    /* We initialize this vector. */

    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
    {
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        if (fe.has_support_on_face (i, f))
        {
          m_fe_sup_on_face[f].push_back (i);
          //std::cout << "m_fe_sup_on_face[f].push_back (i) = " << f << "  "
          //    << i << std::endl; /* bugs*/
        }
      }
    }

    /* Initialize the memory for the vectors. */
    init ();
  }

  template <int dim>
  void
  BoundaryValues<dim>::init ()
  {
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const unsigned int n_angles = mp_ord->get_n_angles ();

    /* Fill the containers with the incoming faces and owners. */
    m_i_faces.resize (n_angles);
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      set_incoming_faces (i);
    }

    /* Fill the containers with the outgoing faces and owners. */
    m_o_faces.resize (n_angles);
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      set_outgoing_faces (i);
    }

    /* This is a check. If this is not fulfilled, something is wrong. */
    {
      /* Count total number of incoming faces per energy group. */
      unsigned int i_faces_per_group = 0;
      unsigned int o_faces_per_group = 0;
      for (unsigned int i = 0; i < n_angles; ++i)
      {
        i_faces_per_group += m_i_faces[i].size ();
        o_faces_per_group += m_o_faces[i].size ();
      }
      /* Incoming and outgoing should be the same number. */
      Assert(i_faces_per_group == o_faces_per_group,
          ExcMessage("Incoming faces different from outgoing faces."));

      /* After the check, we liberate the space used by m_o_faces. */
      m_o_faces.clear ();
    }

    /* Number of dofs per direction for the boundary conditions data. */
    const unsigned int dofs_per_face = std::pow(fe.degree + 1, dim-1);
    m_dofs_in_bc.resize (n_angles);
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      m_dofs_in_bc[i] = m_i_faces[i].size () * dofs_per_face;
    }

    /* Set the global_to_bc_data vector. */
    m_g2s_dofs.resize (n_angles);
    for (unsigned int i = 0; i < n_angles; ++i)
    {
      set_global_to_bc_data (i);
    }
  }

  template <int dim>
  void
  BoundaryValues<dim>::set_incoming_faces (const unsigned int i_ord)
  {
    std::multimap<c_iter, unsigned int> & b_faces = m_i_faces[i_ord];

    Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const MappingQ1<dim> & mapping (mp_state.get_mapping ());

    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_quadrature_points
        | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mapping, fe, face_quadrature,
        face_update_flags);

    /* clear the container before adding new elements with insert. */
    b_faces.clear ();

    typename DoFHandler<dim, dim>::active_cell_iterator cell, endc;
    cell = mp_state.get_dof_handler ().begin_active ();
    endc = mp_state.get_dof_handler ().end ();
    for (; cell != endc; ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        /* face is the element we are interested. */
        f_iter face = cell->face (f);

        /* if face is not in a boundary we are not interested. */
        if (not face->at_boundary ())
        {
          continue;
        }

        /* If boundary face, but is not albedo, not interested. */
        /* unsigned int fbi = face->boundary_id ();
         if (mp_state.mp_data.mp_geom.core.bcs[fbi] != 2)
         {
         continue;
         } */

        /* Get normal outward vector to the face. */
        fe_values_face.reinit (cell, f);

        /* if we are parallel to the boundary we are not interested. */
        const Tensor<1, dim> &n0 = fe_values_face.normal_vector (0);
        if (std::abs (beta * n0) < 1.e-10)
        {
          continue;
        }

        /* if it is an outgoing boundary we are not interested. */
        if (beta * n0 > 0)
        {
          continue;
        }

        /* if passing all tests, the face is incoming or albedo boundary. */
        b_faces.insert (std::pair<c_iter, unsigned int> (cell, f));
      }
    }
  }

  template <int dim>
  void
  BoundaryValues<dim>::set_outgoing_faces (const unsigned int i_ord)
  {
    std::multimap<c_iter, unsigned int> & b_faces = m_o_faces[i_ord];

    Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_ord));
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    const MappingQ1<dim> & mapping (mp_state.get_mapping ());

    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_quadrature_points
        | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mapping, fe, face_quadrature,
        face_update_flags);

    // clear the container before adding new elements with insert
    b_faces.clear ();

    typename DoFHandler<dim, dim>::active_cell_iterator cell, endc;
    cell = mp_state.get_dof_handler ().begin_active ();
    endc = mp_state.get_dof_handler ().end ();
    for (; cell != endc; ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        /* Face is the element we are interested in. */
        f_iter face = cell->face (f);

        /* If face is not in a boundary we are not interested. */
        if (not face->at_boundary ())
          continue;

        /* If boundary face, but is not incoming (or albedo), not interested. */
        /* unsigned int fbi = face->boundary_id ();
         if (mp_state.mp_data.mp_geom.core.bcs[fbi] != 2)
         {
         continue;
         } */

        /* Get normal outward vector to the face. */
        fe_values_face.reinit (cell, f);

        const Tensor<1, dim> &n0 = fe_values_face.normal_vector (0);
        /* if we are parallel to the boundary we are not interested. */
        if (std::abs (beta * n0) < 1.e-10)
        {
          continue;
        }

        /* If it is an incoming boundary we are not interested. */
        if (beta * n0 < 0)
        {
          continue;
        }

        // after passing all the previous tests, the face should be an
        // incoming or albedo boundary for the present direction, so we add
        // it to the container
        b_faces.insert (std::pair<c_iter, unsigned int> (cell, f));
      }
    }
  }

  template <int dim>
  void
  BoundaryValues<dim>::set_global_to_bc_data (const unsigned int i_ord)
  {

    std::map<unsigned int, unsigned int> & map_dofs = m_g2s_dofs[i_ord];
    std::multimap<c_iter, unsigned int> & b_faces = m_i_faces[i_ord];

    // clear the container before adding new elements with push_back
    map_dofs.clear ();

    /* Temporary vector with the local dofs indices. */
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    std::vector<types::global_dof_index> l_dof_ind (fe.dofs_per_cell);

    unsigned int bc_ind = 0; /* A different index for each dofs in the bc. */
    typename std::multimap<c_iter, unsigned int>::iterator it;
    for (it = b_faces.begin (); it != b_faces.end (); ++it)
    {
      c_iter cell = it->first;
      unsigned int f = it->second;
      f_iter face = cell->face (f);
      cell->get_dof_indices (l_dof_ind);

      for (unsigned int i = 0; i < m_fe_sup_on_face[f].size (); ++i, ++bc_ind)
      {
        const unsigned int ii = m_fe_sup_on_face[f][i];
        /* map_dofs.insert (
         *  std::pair<unsigned int, unsigned int> (l_dof_ind[ii], bc_ind)); */
        /* make_pair implicitly deduce the type of the arguments. */
        map_dofs.insert (std::make_pair (l_dof_ind[ii], bc_ind));
      }
    }
  }

  template <int dim>
  void
  BoundaryValues<dim>::apply_reflective (Vector<double> & src_m,
                                         const BlockVector<double> & psi_s,
                                         const unsigned int /* g */,
                                         const unsigned int i_dir,
                                         const unsigned int f_bc_ind)
  {
    Point<dim> beta = vector_to_point<dim> (mp_ord->get_q (i_dir));

    // we get the finite elements
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());

    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
                                          | update_quadrature_points
                                          | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mp_state.get_mapping (), fe,
        face_quadrature, face_update_flags);

    const unsigned int dofs_per_face = std::pow(fe.degree + 1, dim-1);
    Vector<double> cell_src (dofs_per_face), cell_dst (dofs_per_face);
    std::vector<unsigned int> l_dof_ind (fe.dofs_per_cell);

    typename std::multimap<c_iter, unsigned int>::iterator it;
    for (it = m_i_faces[i_dir].begin (); it != m_i_faces[i_dir].end (); ++it)
    {
      c_iter cell = it->first;
      unsigned int f = it->second;
      f_iter face = cell->face (f);

      fe_values_face.reinit (cell, f);
      const std::vector<double> &JxW = fe_values_face.get_JxW_values ();

      /* Get normal outward vector to the face. */
      const Tensor<1, dim> &n0 = fe_values_face.normal_vector (0);
      double alb_bv = 1.0;
      //-----------------------------------------------------------------------
      /**
       @todo use extract_bv to set this albedos
       @code{.cpp}
       if (mp_state.mp_data.mp_settings.get_bcs_from_file ())
       {
       // here I should look for the albedo coefficient!!!
       unsigned int ass_i = mp_state.mp_geom.get_cells_to_assembly (
       cell->user_index ());
       unsigned int pin_i = mp_state.mp_geom.get_cells_to_pin (
       cell->user_index ());
       unsigned int f = -1;
       for (unsigned int ref_f = 0; ref_f < GeometryInfo<dim>::faces_per_cell;
       ++ref_f)
       {
       if (std::abs (p_normal<dim> (ref_f) * n0 - 1.) < 1.e-6)
       {
       f = ref_f;
       }
       }
       unsigned int loc_pin =
       mp_state.m_assembly_boundary_pins[ass_i][f][pin_i];
       alb_bv = mp_state.m_albedo_bv[ass_i][f][g][i_ord][loc_pin];
       }
       @endcode
       */
      //-----------------------------------------------------------------------
      cell->get_dof_indices (l_dof_ind); /* local dof indices*/
      cell_dst = 0;

      // The faces are straight lines, so the normals are independent of
      // the quadrature point
      const double beta_n = beta * n0;
      for (unsigned int q = 0; q < fe_values_face.n_quadrature_points; ++q)
      {
        double cell_src_q = 0;
        for (unsigned int j = 0; j < m_fe_sup_on_face[f].size (); ++j)
        {
          const unsigned int jj = m_fe_sup_on_face[f][j];
          cell_src_q += psi_s.block (f_bc_ind + i_dir) (
                            m_g2s_dofs[i_dir][l_dof_ind[jj]])
                        * fe_values_face.shape_value (jj, q);
        }
        for (unsigned int i = 0; i < m_fe_sup_on_face[f].size (); ++i)
        {
          const unsigned int ii = m_fe_sup_on_face[f][i];
          cell_dst (i) -= beta_n * cell_src_q * alb_bv
                          * fe_values_face.shape_value (ii, q) * JxW[q];
        }
      }
      for (unsigned int i = 0; i < m_fe_sup_on_face[f].size (); ++i)
      {
        const unsigned int ii = m_fe_sup_on_face[f][i];
        src_m[l_dof_ind[ii]] += cell_dst (i);
      }
    }
  }

  template <int dim>
  void
  BoundaryValues<dim>::get_reflective (BlockVector<double> & psi_s,
                                       const Vector<double> & psi_m,
                                       const unsigned int o_dir,
                                       const unsigned int f_bc_ind)
  {

    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    /* the dofs_per_face is not well defined in 1d, so we set it up manually. */
    const unsigned int dofs_per_face = std::pow(fe.degree + 1, dim-1);

    std::vector<types::global_dof_index> l_dof_ind (fe.dofs_per_cell);
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
    {
      /* Identify incoming dir. i_dir associated with o_dir at each face. */
      unsigned int i_dir = mp_ord->out2in[f][o_dir];

      /* The ougoing flux is not inserted back by reflection. */
      if (i_dir == (unsigned int) -1)
      {
        continue;
      }

      /* This face, f, does not have reflecting bcs. */
      if (mp_state.mp_data.mp_geom.m_core.m_bcs[f] == 2)
      {
        typename std::multimap<c_iter, unsigned int>::iterator it;
        for (it = m_i_faces[i_dir].begin (); it != m_i_faces[i_dir].end (); ++it)
        {
          c_iter cell = it->first;
          unsigned int face_no = it->second;
          f_iter face = cell->face (face_no);
          /* @todo I am not sure about this...? */
          if (face->boundary_id () != f)
          {
            continue;
          }
          cell->get_dof_indices (l_dof_ind);
          for (unsigned int i = 0; i < dofs_per_face; ++i)
          {
            const unsigned int ii = m_fe_sup_on_face[face_no][i];
            psi_s.block (f_bc_ind + i_dir)[m_g2s_dofs[i_dir][l_dof_ind[ii]]] =
                psi_m[l_dof_ind[ii]];
          }
        }
      }
      else
      {
        continue;
        // If it is not a reflective boundary it should be "zero incoming flux"
        // until I add the option of "given incoming flux"
        //psi_s.block (f_bc_ind + i_dir) = 0.;
      }
    }
  }

  template class BoundaryValues<1> ;
  template class BoundaryValues<2> ;
  template class BoundaryValues<3> ;

} // end of namespace Forest
