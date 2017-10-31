/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/extract_bv.cc
 * @brief  Implementation of class template ExtractBV
 */

#include "neutronics/extract_bv.h"

#include "geometry/prob_geom.h"
#include "input/input.h"
#include "input/input_settings.h"
#include "input/input_mat.h"
#include "utils/forest_utils_dealii.h"
#include "geometry/path_to_gridgenerator_dealii.h"
#include "neutronics/downstream_new.h"
#include "neutronics/state.h"
#include "angle/moment_indexer.h"
#include "angle/spherical_harmonics.h"
#include "angle/quadrature_base.h"
#include "angle/quadrature_factory.h"
#include "angle/moments_to_directions.h"

#include "deal.II/base/config.h"        // for DEAL_II_VERSION_GTE
#if DEAL_II_VERSION_GTE(8,4,0)
#include "deal.II/base/tensor.h"        // for Tensor
#endif
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/point.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/grid/tria.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <cmath>

namespace Forest
{

  template <int dim>
  ExtractBV<dim>::ExtractBV (const State<dim> &state)
  : m_verbose (false),
  mp_state (state),
  m_ord (QuadratureFactory<dim>::build (
          mp_state.mp_data.mp_settings.get_sn_order (),
          mp_state.mp_data.mp_settings.get_quad_name ())),
  mp_geom (mp_state.mp_geom),
  m_dof_handler (mp_state.get_dof_handler ()),
  m_n_groups (mp_state.get_n_groups ()),
  m_n_angles (m_ord->get_n_angles ())
  {
    m_bcs_order.resize (1, 0);
    this->initialize ();
  }

  template <int dim>
  void
  ExtractBV<dim>::initialize ()
  {

    // This three functions must be called in this order, so I must lock it.
    this->set_pin_bc ();// Set the container for the pin boundaries: (cell, f)
    this->set_pin_bw ();// Set the weights for each pin face
    this->set_pin_bv ();// Extract the boundary flux

    this->set_ass_bc ();// Set the container for the ass boundaries: (ass, f)
    this->set_ass2pin ();// Set pointers (assembly,face) to (pin, counter)
    this->set_ass_bw ();// Set the weights for each ass face
    this->set_ass_bv ();// Extract the boundary flux

    // Calculate the values for the partial currents and fluxes over boundaries
    this->calculate_current_pin ();
    this->calculate_current_ass ();

    // This is the new information for the boundary values (at least for 1d)
    this->calculate_moments_pin ();
    this->calculate_moments_ass ();

    // We print some information if required
    if (m_verbose)
    this->print_ass_bc ();

    // Read values from file?
    /*
     if (mp_state.mp_data.mp_settings.get_bcs_from_file ())
     this->read_ass_bv ();
     */

  }

  template <int dim>
  void
  ExtractBV<dim>::set_pin_bc ()
  {
    // We use some constants
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // Alias for the map
    //typename std::map<c_iter, unsigned int> map_cf;
    using map_cf = std::map<c_iter, unsigned int>;

    // We resize the container with
    m_pin_bc.resize (n_pins, std::vector<map_cf> (n_faces));

    // We use the finite element to obtain the normal vector
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_quadrature_points
    | update_normal_vectors;
    FEFaceValues<dim> fe_vface (fe, face_quadrature, face_update_flags);

    // We define the iterators an run over all the cells
    c_iter cell_it = m_dof_handler.begin_active (), endc = m_dof_handler.end ();
    for (; cell_it != endc; ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = mp_geom.get_cell2pin (cell);
      // Now we run over all the faces
      for (unsigned int f = 0; f < n_faces; ++f)
      {
        if (cell_it->face(f)->at_boundary())
        {
          fe_vface.reinit (cell_it, f);
          const Tensor<1,dim> &normal0 = fe_vface.normal_vector (0);
          const unsigned int face = normal2face (normal0);
          m_pin_bc[pin][face].insert (std::make_pair (cell_it, f));
          continue;
        }
        unsigned int cell_neigh;
        if (cell_it->neighbor(f)->has_children())
        {
          cell_neigh = cell_it->neighbor(f)->child_index(0);
        }
        else
        {
          cell_neigh = cell_it->neighbor(f)->user_index();
        }
        const unsigned int pin_neigh = mp_geom.get_cell2pin (cell_neigh);
        // If the neighbor cell belong to a different pin, we want its face
        if (pin != pin_neigh)
        {
          fe_vface.reinit (cell_it, f);
#if DEAL_II_VERSION_GTE(8,4,0)
          const Tensor<1,dim> &normal0 = fe_vface.normal_vector (0);
#else
          const Point<dim> &normal0 = fe_vface.normal_vector (0);
#endif
          const unsigned int face = normal2face (normal0);
          m_pin_bc[pin][face].insert (std::make_pair (cell_it, f));
        }
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_ass_bc ()
  {
    // We talk
    if (m_verbose)
    std::cout << " set_ass_boundary " << std::endl;

    // We use some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We resize the container with
    m_ass_bc.resize (n_ass,
        std::vector<std::map<c_iter, unsigned int> > (n_faces));

    // We use the finite element to obtain the normal vector
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_quadrature_points
    | update_normal_vectors;
    FEFaceValues<dim> fe_vface (fe, face_quadrature, face_update_flags);

    // We define the iterators an run over all the cells
    c_iter cell_it = m_dof_handler.begin_active (), endc = m_dof_handler.end ();
    for (; cell_it != endc; ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = mp_geom.get_cell2pin (cell);
      const unsigned int ass = mp_geom.get_pin2assembly (pin);
      // Now we run over all the faces
      for (unsigned int f = 0; f < n_faces; ++f)
      {
        if (cell_it->face(f)->at_boundary())
        {
          fe_vface.reinit (cell_it, f);
#if DEAL_II_VERSION_GTE(8,4,0)
          const Tensor<1,dim> &normal0 = fe_vface.normal_vector (0);
#else
          const Point<dim> &normal0 = fe_vface.normal_vector (0);
#endif
          const unsigned int face = normal2face (normal0);
          m_ass_bc[ass][face].insert (std::make_pair (cell_it, f));
          continue;
        }
        unsigned int cell_neigh;
        if (cell_it->neighbor(f)->has_children())
        {
          cell_neigh = cell_it->neighbor(f)->child_index(0);
        }
        else
        {
          cell_neigh = cell_it->neighbor(f)->user_index();
        }
        const unsigned int ass_neigh = mp_geom.get_cell2assembly (cell_neigh);
        // If the neighbor cell belong to a different assembly, we want its face
        if (ass != ass_neigh)
        {
          fe_vface.reinit (cell_it, f);
#if DEAL_II_VERSION_GTE(8,4,0)
          const Tensor<1,dim> &normal0 = fe_vface.normal_vector (0);
#else
          const Point<dim> &normal0 = fe_vface.normal_vector (0);
#endif
          const unsigned int face = normal2face (normal0);
          m_ass_bc[ass][face].insert (std::make_pair (cell_it, f));
        }
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_ass2pin ()
  {
    // We talk
    if (m_verbose)
    std::cout << " set_assembly_boundary_pins " << std::endl;

    // We use some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We resize the container with
    m_ass2pin.resize (n_ass,
        std::vector<std::set<unsigned int> > (n_faces));

    // We run over all the assemblies and faces of each assembly
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      // We define the iterators an run over all the cells
      typename std::map<c_iter, unsigned int>::iterator it, end;
      it = m_ass_bc[ass][f].begin ();
      end = m_ass_bc[ass][f].end ();
      for (; it != end; ++it)
      {
        // We extract the id of the pin the cell belongs to and store it
        const unsigned int cell = it->first->user_index ();
        const unsigned int pin = mp_geom.get_cell2pin (cell);
        m_ass2pin[ass][f].insert (pin);
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::print_ass_bc ()
  {
    // We talk
    if (m_verbose)
    std::cout << " print_ass_boundary " << std::endl;

    // We use some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We run over all the assemblies and faces of each assembly
    if (m_verbose)
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      std::cout << "ass , f = " << ass << " , " << f << std::endl;
      typename std::map<c_iter, unsigned int>::iterator it, end;
      it = m_ass_bc[ass][f].begin ();
      end = m_ass_bc[ass][f].end ();
      for (; it != end; ++it)
      {
        const unsigned int cell = it->first->user_index ();
        const unsigned int face = it->second;
        std::cout << "  cell,face = " << cell << "," << face << std::endl;
      }
    }
    if (m_verbose)
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      std::cout << "m_ass_bc[ass][f].size() = " << m_ass_bc[ass][f].size ()
      << std::endl;
      typename std::map<c_iter, unsigned int>::iterator it, end;
      it = m_ass_bc[ass][f].begin ();
      end = m_ass_bc[ass][f].end ();
      for (; it != end; ++it)
      {
        const unsigned int cell = it->first->user_index ();
        const unsigned int face = it->second;
        std::cout << "  cell,face = " << cell << "," << face << std::endl;
      }
    }
    if (m_verbose)
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      std::cout << "m_ass2pin[ass][f].size() = "
      << m_ass2pin[ass][f].size () << std::endl;
      typename std::set<unsigned int>::iterator it, end;
      it = m_ass2pin[ass][f].begin ();
      end = m_ass2pin[ass][f].end ();
      for (; it != end; ++it)
      {
        const unsigned int pin = *it;
        std::cout << "  pin = " << pin << std::endl;
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_pin_bw ()
  {
    // We need some constants
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We resize the containers
    m_pin_bw.resize (n_pins, std::vector<double> (n_faces, 0.));

    // Finite element things
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
    | update_quadrature_points
    | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mp_state.get_mapping (), fe,
        face_quadrature, face_update_flags);
    const unsigned int n_q = fe_values_face.n_quadrature_points;

    // Run over the pin faces to extract average values per group and angle
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      // In order to store average spatial values, we calculate the measure
      // of the pin face to divide the integral of the flux by this value
      m_pin_bw[pin][f] = 0.;

      // We run over the cells composing this pin face and integrate
      typename map_cf::iterator it, end;
      it = m_pin_bc[pin][f].begin ();
      end = m_pin_bc[pin][f].end ();
      for (; it != end; ++it)
      {
        c_iter cell_it = it->first;
        const unsigned int f2 = it->second;

        // Update fem face information
        fe_values_face.reinit (cell_it, f2);
        const std::vector<double> &JxW = fe_values_face.get_JxW_values ();

        // Cumulative contribution of the cell face length
        double cell_f_measure = 0;
        for (unsigned int q = 0; q < n_q; ++q)
        cell_f_measure += JxW[q];
        m_pin_bw[pin][f] += cell_f_measure;
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_pin_bv ()
  {
    // We need some constants
    //const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_groups = mp_state.get_n_groups ();

    // Resize the container
    m_pin_bv.resize (n_pins,
        stdvd3 (n_faces, stdvd2 (n_groups, stdvd1 (m_n_angles, 0.))));

    // Approximate psi_n and then extract the boundary values with fill_bv_g_n
    for (unsigned int g = 0; g < n_groups; ++g)
    for (unsigned int n = 0; n < m_n_angles; ++n)
    {
      // pointer to the angular flux for group g and direction n
      const Vector<double> & psi_g_n = mp_state.m_psi.m_data[g].block(n);
      // get the bv for this group and angle
      this->set_pin_bv_g_n (psi_g_n, g, n);
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_pin_bv_g_n (const Vector<double> & psi_n,
      const unsigned int g,
      const unsigned int n)
  {
    // We need some constants
    const unsigned int n_pins = mp_geom.get_n_pins ();

    // Finite element things
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
    | update_quadrature_points
    | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mp_state.get_mapping (), fe,
        face_quadrature, face_update_flags);
    const unsigned int n_q = fe_values_face.n_quadrature_points;
    std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);

    // Run over the pin faces to extract average values per group and angle
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
    {
      // We run over the cells composing this pin face and integrate
      typename map_cf::iterator it, end;
      it = m_pin_bc[pin][f].begin ();
      end = m_pin_bc[pin][f].end ();
      for (; it != end; ++it)
      {
        c_iter cell_it = it->first;
        const unsigned int f2 = it->second;

        // Update fem face information
        fe_values_face.reinit (cell_it, f2);
        const std::vector<double> &JxW = fe_values_face.get_JxW_values ();
        cell_it->get_dof_indices (local_dof_indices);

        // Cumulative contribution of the flux per cell face
        double value = 0;
        for (unsigned int q = 0; q < n_q; ++q)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
        value += psi_n (local_dof_indices[j])
        * fe_values_face.shape_value (j, q) * JxW[q];

        //[pin][face][group_g][angle]
        m_pin_bv[pin][f][g][n] += value;
      }
      // divide by the length of the interval to obtain the normalized value
      m_pin_bv[pin][f][g][n] = m_pin_bv[pin][f][g][n] / m_pin_bw[pin][f];
      Assert(numbers::is_finite(m_pin_bv[pin][f][g][n]),
             ExcMessage("Not finite number found"));
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_ass_bw ()
  {
    // We need some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We resize the containers
    m_ass_bw.resize (n_ass, std::vector<double> (n_faces, 0.));

    // Finite element things
    const FiniteElement<dim, dim> & fe (mp_state.get_fe ());
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    const UpdateFlags face_update_flags = update_values | update_JxW_values
    | update_quadrature_points
    | update_normal_vectors;
    FEFaceValues<dim> fe_values_face (mp_state.get_mapping (), fe,
        face_quadrature, face_update_flags);
    const unsigned int n_q = fe_values_face.n_quadrature_points;

    // Run over the pin faces to extract average values per group and angle
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      // In order to store average spatial values, we calculate the measure
      // of the ass face to divide the integral of the flux by this value
      m_ass_bw[ass][f] = 0.;

      // We run over the cells composing this ass face and integrate
      typename map_cf::iterator it, end;
      it = m_ass_bc[ass][f].begin ();
      end = m_ass_bc[ass][f].end ();
      for (; it != end; ++it)
      {
        c_iter cell_it = it->first;
        const unsigned int f2 = it->second;

        // Update fem face information
        fe_values_face.reinit (cell_it, f2);
        const std::vector<double> &JxW = fe_values_face.get_JxW_values ();

        // Cumulative contribution of the cell face length
        double cell_f_measure = 0;
        for (unsigned int q = 0; q < n_q; ++q)
        cell_f_measure += JxW[q];

        m_ass_bw[ass][f] += cell_f_measure;
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::set_ass_bv ()
  {
    // We need some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_groups = mp_state.get_n_groups ();

    // Resize the container
    m_ass_bv.resize (n_ass,
        stdvd3 (n_faces, stdvd2 (n_groups, stdvd1 (m_n_angles, 0.))));

    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      typename std::set<unsigned int>::iterator it, end;
      it = m_ass2pin[ass][f].begin ();
      end = m_ass2pin[ass][f].end ();
      for (; it != end; ++it)
      {
        const unsigned int pin = *it;
        const double weight = m_pin_bw[pin][f] / m_ass_bw[ass][f];
        for (unsigned int g = 0; g < n_groups; ++g)
        for (unsigned int n = 0; n < m_n_angles; ++n)
        m_ass_bv[ass][f][g][n] += m_pin_bv[pin][f][g][n] * weight;
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::calculate_current_pin ()
  {
    // We need some constants
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We need some constants for the quadrature
    const std::vector<double> & w = m_ord->get_w ();
    const std::vector<std::vector<double> > & q = m_ord->get_q ();
    const double eps = 1.e-10;

    // Resize the containers and initialize to zero
    m_phi0_pin.resize (n_pins, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jp_o_pin.resize (n_pins, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jp_i_pin.resize (n_pins, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jnet_pin.resize (n_pins,
        stdvd3 (n_faces, stdvd2 (n_groups, stdvd1 (dim, 0.))));
    m_j_pin.resize (n_pins, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_phi2_pin.resize (n_pins, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));

    for (unsigned int pin = 0; pin < n_pins; ++pin)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      Point<dim> p_normal (p_normal<dim> (f));
      const double normal_norm = std::sqrt (p_normal.square ());
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        for (unsigned int n = 0; n < m_n_angles; ++n)
        {
          // Generate the scalar flux phi0
          m_phi0_pin[pin][f][g] += w[n] * m_pin_bv[pin][f][g][n];
          // Generate the net current vector
          for (unsigned int d = 0; d < dim; ++d)
          m_jnet_pin[pin][f][g][d] += w[n] * q[n][d] * m_pin_bv[pin][f][g][n];
          // Generate the outgoing and incoming partial current J+ J-
          const double q_dot_n = vector_to_point<dim> (q[n]) * p_normal;
          if (q_dot_n > eps)
          m_jp_o_pin[pin][f][g] += w[n] * (q_dot_n) * m_pin_bv[pin][f][g][n];
          else
          m_jp_i_pin[pin][f][g] += w[n] * (-q_dot_n) * m_pin_bv[pin][f][g][n];
        }
        m_j_pin[pin][f][g] = m_jp_o_pin[pin][f][g] - m_jp_i_pin[pin][f][g];

        // cosinus between the net current and the outgoing normal
        Point<dim> p_jnet (vector_to_point<dim> (m_jnet_pin[pin][f][g]));
        const double jnet_dot_norm = p_jnet * p_normal;
        const double jnet_norm = std::sqrt (p_jnet.square ());
        double costheta = 0.0;
        // if theta == 0, then costheta^2 = 1, then p2(costheta) = 1
        if (jnet_norm < eps)
          costheta = 1;
        else
          costheta = jnet_dot_norm / (jnet_norm * normal_norm);
        // second order legendre polynomial 1/2*(3*x^2-1)
        const double p2costheta = 1. / 2 * (3 * std::pow (costheta, 2) - 1.);
        m_phi2_pin[pin][f][g] =
        3. / (5. * p2costheta) * (m_jp_o_pin[pin][f][g] + m_jp_i_pin[pin][f][g]
            - m_phi0_pin[pin][f][g] / 2.);

        /*
         std::cout << "m_jp_o[pin][f][g] = " << std::setprecision(10) << m_jp_o[pin][f][g] << std::endl;
         std::cout << "m_jp_i[pin][f][g] = " << std::setprecision(10) << m_jp_i[pin][f][g] << std::endl;
         std::cout << "p_jnet = " << p_jnet << std::endl;
         std::cout << "jnet_dot_norm = " << std::setprecision(10) << jnet_dot_norm << std::endl;
         std::cout << "jnet_norm = " << std::setprecision(10) << jnet_norm << std::endl;
         std::cout << "costheta = " << costheta << std::endl;
         std::cout << "m_phi2[pin][f][g] = " << m_phi2[pin][f][g] << std::endl;
         std::cout << std::endl;*/
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::calculate_current_ass ()
  {
    // We need some constants
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // We need some constants for the quadrature
    const std::vector<double> & w = m_ord->get_w ();
    const std::vector<std::vector<double> > & q = m_ord->get_q ();
    const double eps = 1.e-10;

    // Resize the containers and initialize to zero
    m_phi0_ass.resize (n_ass, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jp_o_ass.resize (n_ass, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jp_i_ass.resize (n_ass, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_jnet_ass.resize (n_ass,
        stdvd3 (n_faces, stdvd2 (n_groups, stdvd1 (dim, 0.))));
    m_j_ass.resize (n_ass, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));
    m_phi2_ass.resize (n_ass, stdvd2 (n_faces, stdvd1 (n_groups, 0.)));

    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      Point<dim> p_normal (p_normal<dim> (f));
      const double normal_norm = std::sqrt (p_normal.square ());
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        for (unsigned int n = 0; n < m_n_angles; ++n)
        {
          // Generate the scalar flux phi0
          m_phi0_ass[ass][f][g] += w[n] * m_ass_bv[ass][f][g][n];
          // Generate the net current vector
          for (unsigned int d = 0; d < dim; ++d)
          m_jnet_ass[ass][f][g][d] += w[n] * q[n][d] * m_ass_bv[ass][f][g][n];
          // Generate the outgoing and incoming partial current J+ J-
          const double q_dot_n = vector_to_point<dim> (q[n]) * p_normal;
          if (q_dot_n > eps)
          m_jp_o_ass[ass][f][g] += w[n] * (q_dot_n) * m_ass_bv[ass][f][g][n];
          else
          m_jp_i_ass[ass][f][g] += w[n] * (-q_dot_n) * m_ass_bv[ass][f][g][n];
        }
        m_j_ass[ass][f][g] = m_jp_o_ass[ass][f][g] - m_jp_i_ass[ass][f][g];

        // cosinus between the net current and the outgoing normal
        Point<dim> p_jnet (vector_to_point<dim> (m_jnet_ass[ass][f][g]));
        const double jnet_dot_norm = p_jnet * p_normal;
        const double jnet_norm = std::sqrt (p_jnet.square ());
        double costheta = 0.0;
        // if theta == 0, then costheta^2 = 1, then p2(costheta) = 1
        if (jnet_norm < eps)
          costheta = 1;
        else
          costheta = jnet_dot_norm / (jnet_norm * normal_norm);
        // second order legendre polynomial 1/2*(3*x^2-1)
        const double p2costheta = 1. / 2 * (3 * std::pow (costheta, 2) - 1.);
        m_phi2_ass[ass][f][g] =
        3. / (5. * p2costheta) * (m_jp_o_ass[ass][f][g] + m_jp_i_ass[ass][f][g]
            - m_phi0_ass[ass][f][g] / 2.);
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::calculate_moments_pin ()
  {
    // We need some constants
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    const unsigned int leg_mom = 5;

    // We need some constants for the quadrature
    const std::vector<double> & w = m_ord->get_w ();
    const std::vector<std::vector<double> > & q = m_ord->get_q ();

    MomentIndexer<dim> ind (leg_mom);
    const unsigned int n_moments = ind.get_number_of_moments ();

    // Resize the containers and initialize to zero
    m_pin_bvm.resize (n_pins, stdvd3(n_moments, stdvd2 (n_faces, stdvd1 (n_groups,0.))));

    // run over the pins, faces and groups to perform angular coarsening
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    for (unsigned int f = 0; f < n_faces; ++f)
    for (unsigned int g = 0; g < m_n_groups; ++g)
    {
      // short for the discrete ordinates representation of the flux
      std::vector<double> & flux_directions = m_pin_bv[pin][f][g];

      // we project the angular flux into the spherical harmonics moments
      for (unsigned int m = 0; m < n_moments; ++m)
      {
        std::vector<int> lm =
        { ind.get_l (m), ind.get_m (m)};
        for (unsigned int n = 0; n < m_n_angles; ++n)
        {
          const double y_lm = SphericalHarmonics<dim>::Y_lm (lm, q[n]);
          m_pin_bvm[pin][m][f][g] += w[n] * y_lm * flux_directions[n];
        }
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::calculate_moments_ass ()
  {
    // We need some constants
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    const unsigned int leg_mom = 5;

    // We need some constants for the quadrature
    const std::vector<double> & w = m_ord->get_w ();
    const std::vector<std::vector<double> > & q = m_ord->get_q ();

    MomentIndexer<dim> ind (leg_mom);
    const unsigned int n_moments = ind.get_number_of_moments ();

    // Resize the containers and initialize to zero
    m_ass_bvm.resize (n_ass, stdvd3(n_moments, stdvd2 (n_faces, stdvd1 (n_groups,0.))));

    // run over the pins, faces and groups to perform angular coarsening
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    for (unsigned int f = 0; f < n_faces; ++f)
    for (unsigned int g = 0; g < m_n_groups; ++g)
    {
      // short for the discrete ordinates representation of the flux
      std::vector<double> & flux_directions = m_ass_bv[ass][f][g];

      // we project the angular flux into the spherical harmonics moments
      for (unsigned int m = 0; m < n_moments; ++m)
      {
        std::vector<int> lm =
        { ind.get_l (m), ind.get_m (m)};
        for (unsigned int n = 0; n < m_n_angles; ++n)
        {
          const double y_lm = SphericalHarmonics<dim>::Y_lm (lm, q[n]);
          m_ass_bvm[ass][m][f][g] += w[n] * y_lm * flux_directions[n];
        }
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::print_ass_bv (const std::string &filename_prefix) const
  {
    // We use some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    for (unsigned int ass = 0; ass < n_ass; ++ass)
    {
      // Create empty property tree object
      using boost::property_tree::ptree;
      ptree pt;
      std::string bin;// bin is a container to save temporal data
      std::string indent ("\t");
      pt.add ("boundary_conditions", "");
      pt.put ("boundary_conditions.<xmlattr>.SN_order", m_ord->get_order ());
      pt.add ("boundary_conditions.keff", mp_state.get_keff ());
      for (unsigned int f = 0; f < n_faces; ++f)
      {
        ptree & node_face = pt.add ("boundary_conditions.face", "");
        node_face.put ("<xmlattr>.id", f);
        for (unsigned int g = 0; g < m_n_groups; ++g)
        {
          ptree & node_g = node_face.add ("group", "");
          node_g.put ("<xmlattr>.id", g);
          vector_to_string (m_ass_bv[ass][f][g], bin, indent, 4);
          node_g.add ("ANGxPOS", bin);
          bin.clear ();
        }
      }
      // Write property tree to XML file (\t is a tabulator)
#if BOOST_VERSION < 105500
      boost::property_tree::xml_writer_settings<char> settings (' ', 2);
#else
      boost::property_tree::xml_writer_settings<std::string> settings (' ', 2);
#endif
      std::string filename = filename_prefix + std::string ("_bc_ass_")
      + num_to_string (ass) + std::string (".xml");
      std::cout << "priting to <" << filename << "> ..." << std::endl;
      write_xml (filename, pt, std::locale (), settings);
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::read_ass_bv ()
  {
    // We talk
    if (m_verbose)
    std::cout << " read_bv " << std::endl;

    // We use some constants
    const unsigned int n_ass = mp_geom.get_n_assemblies ();
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const std::string f_bcs (mp_state.mp_data.mp_settings.get_file_bcs ());

    // We resize the containers
    m_ass_bv.resize (n_ass, stdvd3 (n_faces, stdvd2 (m_n_groups)));
    m_pin_bv.resize (n_pins, stdvd3 (n_faces, stdvd2 (m_n_groups)));

    // We verify that it is only one assembly
    Assert(n_ass == 1, ExcMessage("More than one assembly provided."));

    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    // Load XML file and put its contents in property tree.
    read_xml (f_bcs, pt, boost::property_tree::xml_parser::trim_whitespace);

    std::string bin;// bin is a container to save temporal data
    std::string indent ("\t");

    m_bcs_order[0] = pt.get<unsigned int> (
        "boundary_conditions.<xmlattr>.SN_order");

    // assemblies
    for (const std::pair<std::string, ptree> & face : pt.get_child("boundary_conditions"))
    {
      if (face.first == "face")
      {
        for (const std::pair<std::string, ptree> & group : face.second.get_child(""))
        {
          if (group.first == "group")
          {
            unsigned int f = face.second.get<unsigned int> ("<xmlattr>.id");
            unsigned int g = group.second.get<unsigned int> ("<xmlattr>.id");
            bin = group.second.get<std::string> ("ANGxPOS");
            Forest::string_to_vector (bin, m_ass_bv[0][f][g]);
            bin.clear ();
          }
        }
      }
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::coarsening_spatial ()
  {
    // We talk
    if (m_verbose)
    std::cout << " coarsening_spatial " << std::endl;

    // We use some constants
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // run over the pins, faces and groups to perform angular coarsening
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    for (unsigned int f = 0; f < n_faces; ++f)
    {
      const unsigned int ass = mp_geom.get_pin2assembly (pin);
      for (unsigned int g = 0; g < m_n_groups; ++g)
      for (unsigned int n = 0; n < m_n_angles; ++n)
      m_ass_bv[ass][f][g][n] += m_pin_bv[pin][f][g][n];
    }
  }

  template <int dim>
  void
  ExtractBV<dim>::coarsening_angular ()
  {
    // We talk
    if (m_verbose)
    deallog << " coarsening_angular not implemented" << std::endl;
    /*
    // We need some constants
    const unsigned int n_pins = mp_geom.get_n_pins ();
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

    // variables for the angular coarsening
    std::string q_type = (dim == 1 ? "GaussLegendre" : "LevelSymType2");
    QuadraturePtr<dim> sn_coarse = QuadratureFactory<dim>::build (m_bcs_order,
        q_type);
    std::vector<std::vector<double> > q = sn_coarse->get_q ();
    for (unsigned int n = 0; n < m_n_angles; ++n)
    q[n][2] = std::sqrt (1 - (q[n][0] * q[n][0] + q[n][1] * q[n][1]));
    std::vector<double> w = sn_coarse->get_w ();

    // How many coarsenings in angle?
    const unsigned int ca = mp_state.mp_data.mp_settings.get_coarse_angular ();
    Assert(ca <= sn_coarse->get_n_leg_moments (),
        ExcMessage("Coarsening larger that n_leg_moments provided."));

    // If the coarsening is different from zero, then we project the angular
    // distribution into a spherical harmonics bases of smaller order, and then
    // we project back into the discrete ordinates representation
    if (ca != 0)
    {
      // Moment indexer to run over the spherical harmonics moments
      const unsigned int n_leg_mom_coarse = sn_coarse->get_n_leg_moments ()
      - ca;
      MomentIndexer<dim> ind (n_leg_mom_coarse);
      const unsigned int n_moments = ind.get_number_of_moments ();

      // run over the pins, faces and groups to perform angular coarsening
      for (unsigned int pin = 0; pin < n_pins; ++pin)
      for (unsigned int f = 0; f < n_faces; ++f)
      for (unsigned int g = 0; g < m_n_groups; ++g)
      {
        // short for the discrete ordinates representation of the flux
        std::vector<double> & flux_directions = m_pin_bv[pin][f][g];

        // temporal storage for the angular flux
        std::vector<double> flux_moments (n_moments, 0.);

        // we project the angular flux into the spherical harmonics moments
        for (unsigned int m = 0; m < n_moments; ++m)
        {
          std::vector<int> lm =
          { ind.get_l (m), ind.get_m (m)};
          for (unsigned int n = 0; n < m_n_angles; ++n)
          {
            const double y_lm = SphericalHarmonics<dim>::Y_lm (lm, q[n]);
            flux_moments[m] += w[n] * y_lm * flux_directions[n];
          }
        }

        // clean old values for the angular flux
        for (unsigned int n = 0; n < m_n_angles; ++n)
        flux_directions[n] = 0.;

        // project back to the discrete ordinates representation
        for (unsigned int m = 0; m < n_moments; ++m)
        {
          const std::vector<int> lm =
          { ind.get_l (m), ind.get_m (m)};
          const double scal = (2 * ind.get_l (m) + 1);
          for (unsigned int n = 0; n < m_n_angles; ++n)
          {
            const double y_lm = SphericalHarmonics<dim>::Y_lm (lm, q[n]);
            for (unsigned int pin = 0; pin < n_pins; ++pin)
            flux_directions[n] += y_lm * scal * flux_moments[m];
          }
        }
      }
    }
    */
  }

  template <int dim>
  void
  ExtractBV<dim>::scaling(const double scaling)
  {
    // Some information about the dimensions. The same for pin and ass.
    // double boundary_value = m_pin_bv[pin][f][group_g][angle];
    // double bv = m_jnet[pin][f][group_g][dim];
    // double bv = m_phi0[pin][f][group_g];
    // double bv = m_jp_o[pin][f][group_g];
    // double bv = m_jp_i[pin][f][group_g];
    // double bv = m_j[pin][f][group_g];
    // double bv = m_phi2[pin][f][group_g];

    // Scaling the pin data
    for (unsigned int p = 0; p < m_pin_bv.size(); ++p)
      for (unsigned int f = 0; f < m_pin_bv[p].size(); ++f)
        for (unsigned int g = 0; g < m_pin_bv[p][f].size(); ++g)
        {
          for (unsigned int n = 0; n < m_pin_bv[p][f][g].size(); ++n)
            m_pin_bv[p][f][g][n] *= scaling;
          for (unsigned int d = 0; d < dim; ++d)
            m_jnet_pin[p][f][g][d] *= scaling;

          m_phi0_pin[p][f][g] *= scaling;
          m_jp_o_pin[p][f][g] *= scaling;
          m_jp_i_pin[p][f][g] *= scaling;
          m_j_pin[p][f][g] *= scaling;
          m_phi2_pin[p][f][g] *= scaling;
        }

    // Scaling the assembly data
    for (unsigned int a = 0; a < m_ass_bv.size(); ++a)
      for (unsigned int f = 0; f < m_ass_bv[a].size(); ++f)
        for (unsigned int g = 0; g < m_ass_bv[a][f].size(); ++g)
        {
          for (unsigned int n = 0; n < m_ass_bv[a][f][g].size(); ++n)
            m_ass_bv[a][f][g][n] *= scaling;
          for (unsigned int d = 0; d < dim; ++d)
            m_jnet_ass[a][f][g][d] *= scaling;

          m_phi0_ass[a][f][g] *= scaling;
          m_jp_o_ass[a][f][g] *= scaling;
          m_jp_i_ass[a][f][g] *= scaling;
          m_j_ass[a][f][g] *= scaling;
          m_phi2_ass[a][f][g] *= scaling;
        }
  }

  template <int dim>
  void
  ExtractBV<dim>::scaling_moments(const double scaling)
  {
    // Some information about the dimensions. The same for pin and ass.
    // double boundary_value = m_pin_bv[pin][f][group_g][angle];
    // double bv = m_pin_bvm[pin][n_moments][f][group_g];
    // double bv = m_ass_bvm[ass][n_moments][f][group_g];

    // Scaling the pin data
    for (unsigned int p = 0; p < m_pin_bvm.size(); ++p)
      for (unsigned int m = 0; m < m_pin_bvm[p].size(); ++m)
        for (unsigned int f = 0; f < m_pin_bvm[p][m].size(); ++f)
          for (unsigned int g = 0; g < m_pin_bvm[p][m][f].size(); ++g)
          {
            //deallog << "p,m,f,g = " << p << " " << m << " " << f << " " << g << std::endl;
            m_pin_bvm[p][m][f][g] *= scaling;
          }

    // Scaling the assembly data
    for (unsigned int a = 0; a < m_ass_bvm.size(); ++a)
      for (unsigned int m = 0; m < m_ass_bvm[a].size(); ++m)
        for (unsigned int f = 0; f < m_ass_bvm[a][m].size(); ++f)
          for (unsigned int g = 0; g < m_ass_bvm[a][m][f].size(); ++g)
          {
            //deallog << "a,m,f,g = " << a << " " << m << " " << f << " " << g << std::endl;
            m_ass_bvm[a][m][f][g] *= scaling;
          }
  }

  template class ExtractBV<1> ;
  template class ExtractBV<2> ;
  template class ExtractBV<3> ;

} // end of namespace Forest
