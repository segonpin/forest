/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   geometry/prob_geom.cc
 * @brief  Implementation of class template ProbGeom
 */

#include "geometry/prob_geom.h"
#include "geometry/path_to_gridgenerator_dealii.h"

#include "input/input.h"
#include "geometry/real_core.h"
#include "geometry/real_pin.h"
#include "geometry/real_assembly.h"

#include "utils/forest_utils_logstream.h"
#include "utils/forest_utils_memory.h"
#include "utils/forest_utils_timer.h"

#include "deal.II/base/timer.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/point.h"
#include "deal.II/base/types.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/manifold_lib.h"
#include "deal.II/grid/tria.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/lac/vector.h"
#include "deal.II/numerics/data_out.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <memory>

namespace Forest
{
  using namespace dealii;

  template <>
  void
  ProbGeom<1>::assign_pin_manifolds (Triangulation<1, 1> & /* tria */,
                                     std::vector<RealPin<1> > & /* real_pins */,
                                     unsigned int /* first_manifold */)
  {
    Assert(false, ExcMessage("Not implemented."));
  }

  template <int dim>
  void
  ProbGeom<dim>::assign_pin_manifolds (Triangulation<dim, dim> & tria,
                                       std::vector<RealPin<dim> > & real_pins,
                                       unsigned int first_manifold)
  {
    // After setting the manifold, we relate each cell with its pin
    // in order to apply later different refinements strategies for cells
    // belonging to different type of pins. We do this after each refinement.
    std::vector<unsigned int> cell2pin_tmp;
    cell2pin_tmp = assign_cell_index (tria, real_pins);

    typename Triangulation<dim, dim>::active_cell_iterator cell;
    for (cell = tria.begin_active (); cell != tria.end (); ++cell)
    {
      // We check if we are touching the boundary in the radial direction
      // and we jump this case (no manifold for this cell)
      bool cell_at_boundary = false;
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
      {
        if (cell->face (f)->at_boundary ())
          cell_at_boundary = true;
      }
      if (cell_at_boundary)
        continue;

      const unsigned int pin = cell2pin_tmp[cell->user_index ()];
      if (real_pins[pin].get_type () == PinType::PIN)
      {
        const unsigned int manifold_ind = first_manifold + pin;
        const Point<dim> cell_center = cell->center ();
        // Here 3.0 looks like a magic number, but is the size that we
        // gave to the inner cell at the beginning, and we don't want
        // to refine this cell (or its children).
        const double inner_radius = real_pins[pin].get_radius () / 3.0;
        // we only use x and y coordinates to see if we are inside the pin
        if (real_pins[pin].norm_inf_xy (cell_center) < inner_radius)
        {
          for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
          {
            cell->face (f)->set_all_manifold_ids (manifold_ind);
          }
        }
        else
        {
          //cell->set_all_manifold_ids (types::manifold_id (-1));
        }
      }
      else if (real_pins[pin].get_type () == PinType::PINNEW)
      {
        const unsigned int manifold_ind = first_manifold + pin;
        const Point<dim> cell_center = cell->center ();
        // Here 3.0 looks like a magic number, but is the size that we
        // gave to the inner cell at the beginning, and we don't want
        // to refine this cell (or its children).
        const double inner_radius = real_pins[pin].get_radius () / 3.0;
        // we only use x and y coordinates to see if we are inside the pin
        if (real_pins[pin].norm_inf_xy (cell_center) > inner_radius)
        {
          for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
          {
            if (pin == cell2pin_tmp[cell->neighbor (f)->user_index ()])
              if (cell->material_id () != cell->neighbor (f)->material_id ())
              {
                cell->face (f)->set_all_manifold_ids (manifold_ind);
              }
          }
        }
        else
        {
          cell->set_all_manifold_ids (types::manifold_id (-1));
        }
      }
    }
  }

  template <int dim>
  std::vector<unsigned int>
  ProbGeom<dim>::assign_cell_index (Triangulation<dim, dim> & tria,
                                    std::vector<RealPin<dim> > & real_pins)
  {
    std::vector<unsigned int> cells_to_pin (tria.n_active_cells ());
    unsigned int user_index = 0;
    for (typename Triangulation<dim, dim>::active_cell_iterator cell =
        tria.begin_active (); cell != tria.end (); ++cell)
    {
      for (unsigned int i = 0; i < real_pins.size (); ++i)
      {
        if (real_pins[i].is_inside (cell->center ()))
        {
          cells_to_pin.at (user_index) = i;
          cell->set_user_index (user_index);
          ++user_index;
          // jump to the next cell without finishing the pins loop
          continue;
        }
      }
    }
    return cells_to_pin;
  }

  template <int dim>
  void
  ProbGeom<dim>::refine_assembly (Triangulation<dim, dim> & tria,
                                  std::vector<RealPin<dim> > & real_pins,
                                  const unsigned int n_ref)
  {
    const unsigned int maximum_local_refinements = 3;
    // Here we iterate for the refinements
    for (unsigned int step = 0; step < n_ref; ++step)
    {
      /*if (step >= maximum_local_refinements)
       {
       tria.refine_global (1);
       continue;
       }*/

      // After setting the manifold, we relate each cell with its pin
      // in order to apply later different refinements strategies for cells
      // belonging to different type of pins. We do this after each refinement.
      std::vector<unsigned int> cell2pin_tmp;
      cell2pin_tmp = assign_cell_index (tria, real_pins);

      // Here we look for the faces over the internal boundary, and tag
      // this cells to be refined. (if two faces belong to the internal
      // manifold, will be isotropically refined, otherwise will be refined
      // only in one direction)
      for (typename Triangulation<dim, dim>::active_cell_iterator cell =
          tria.begin_active (); cell != tria.end (); ++cell)
      {
        /* std::cout << "cell->user_index() = " << cell->user_index ();
         Tensor<1, dim> normal = cell->face (4)->center ()
         - cell->face (5)->center ();
         std::cout << "; normal = " << normal << std::endl; */

        // we only use the manifold for the PinType::PIN
        const unsigned int pin = cell2pin_tmp[cell->user_index ()];
        if (real_pins[pin].get_type () == PinType::PIN)
        {
          if (step >= maximum_local_refinements and false)
          {
            const double inner_radius = real_pins[pin].get_radius () * 0.6;
            //if ((real_pins[pin].get_center() - cell->center()).norm() > inner_radius)
            if (real_pins[pin].norm_inf_xy (cell->center ()) > inner_radius)
            {
              //cell->set_refine_flag (RefinementCase<dim>::cut_xy);
              cell->set_refine_flag (RefinementCase<dim>::cut_xy);
            }
          }
          else
          {
            std::vector<unsigned int> direction;
            // we only run over the radial faces (thus, dim = 2 always)
            for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            {
              if (cell->face (f)->manifold_id () != types::manifold_id (-1))
              {
                // Here we use 1-f/2 because:
                // f/2 will be 1 if we want to refine faces f={2,3} (x-axis),
                // then 1-f/2 will be 0, and
                // f/2 will be 0 if we want to refine faces f={0,1} (y-axis),
                // then 1-f/2 will be 1.
                direction.push_back (1 - f / 2);
              }
            }
            if (direction.size () >= 2)
            {
              cell->set_refine_flag (RefinementCase<dim>::cut_xy);
            }
            else if (direction.size () == 1)
            {
              cell->set_refine_flag (
                  RefinementCase<dim>::cut_axis (direction[0]));
            }
            direction.clear ();
          }
        }
        else if (real_pins[pin].get_type () == PinType::PINNEW)
        {
          if (step >= maximum_local_refinements and false)
          {
            const double inner_radius = real_pins[pin].get_radius () / 3.0;
            if (real_pins[pin].norm_inf_xy (cell->center ()) > inner_radius)
            {
              cell->set_refine_flag (RefinementCase<dim>::cut_xy);
            }
          }
          else
          {
            std::vector<unsigned int> direction;
            // we only run over the radial faces (thus, dim = 2 always)
            for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            {
              if (cell->face (f)->manifold_id () != types::manifold_id (-1))
              {
                // Here we use 1-f/2 because:
                // f/2 will be 1 if we want to refine faces f={2,3} (x-axis),
                // then 1-f/2 will be 0, and
                // f/2 will be 0 if we want to refine faces f={0,1} (y-axis),
                // then 1-f/2 will be 1.
                direction.push_back (1 - f / 2);
              }
            }
            if (direction.size () >= 2)
            {
              cell->set_refine_flag (RefinementCase<dim>::cut_xy);
            }
            else if (direction.size () == 1)
            {
              cell->set_refine_flag (
                  RefinementCase<dim>::cut_axis (direction[0]));
            }
            direction.clear ();
          }
        }
        else if (real_pins[pin].get_type () == PinType::PINBOX)
        {
          Point<dim> center = real_pins[pin].get_center ();
          double dx = std::abs ((cell->center () - center)[0]);
          double dy = std::abs ((cell->center () - center)[1]);
          double radius = real_pins[pin].get_radius () * std::sqrt (2.0) / 2.;
          if (dx < radius and dy < radius)
          {
            cell->set_refine_flag (RefinementCase<dim>::cut_xy);
          }
          else
          {
            unsigned int direction = (unsigned int) -1;
            // we look for the face closer to the center, and refine as before
            //double min_dist = std::sqrt(dx*dx + dy*dy);
            for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell;
                ++f, ++f)
            {
              Tensor<1, dim> t0 = (cell->face (f)->center () - center);
              Tensor<1, dim> t1 = (cell->face (f + 1)->center () - center);
              double orthogonality = (t0 * t1) / (t0.norm () * t1.norm ());
              if (std::abs (orthogonality - 1.0) < 1.e-6)
              {
                direction = 1 - f / 2;
                break;
              }
            }
            if (direction != (unsigned int) -1)
            {
              cell->set_refine_flag (RefinementCase<dim>::cut_axis (direction));
            }
          }
        }
        else if (real_pins[pin].get_type () == PinType::BOX)
        {
          if (step >= maximum_local_refinements)
          {
            cell->set_refine_flag (RefinementCase<dim>::cut_xy);
          }
        }
      }
      // there is a bug in 3d, because dealii is adding refinement in the
      // z axis when fixing the refinement of neighbor cells.
      // It can be fixed. I check if something is changed, and in this case
      // I can just run over the cells again and remove the cut_z
      // component just intersecting with the 2d cut_xy flag.
      bool flags_changed = tria.prepare_coarsening_and_refinement ();
      if (flags_changed)
      {
        for (typename Triangulation<dim, dim>::active_cell_iterator cell =
            tria.begin_active (); cell != tria.end (); ++cell)
        {
          RefinementCase<dim> ref_flag = cell->refine_flag_set();
          RefinementCase<dim> ref_xy = RefinementCase<dim>::cut_xy;
          RefinementCase<dim> new_flag = ref_xy & ref_flag;
          cell->set_refine_flag(new_flag);
        }
      }
      // No more than max_loc_refinement anisotropic localized refinements.
      // After that, we do global isotropic refinement.
      tria.execute_coarsening_and_refinement ();
    }
  }

  template <>
  void
  ProbGeom<1>::refine_assembly (Triangulation<1, 1> & tria,
                                std::vector<RealPin<1> > & /* real_pins */,
                                const unsigned int n_ref)
  {
    for (unsigned int step = 0; step < n_ref; ++step)
    {
      tria.refine_global (1);
    }
  }

  /**
   * Set boundary indicators
   */
  template <int dim>
  void
  ProbGeom<dim>::set_boundary_indicators ()
  {
    const double eps = 1.e-10;
    for (typename Triangulation<dim, dim>::active_cell_iterator cell =
        m_tria.begin_active (); cell != m_tria.end (); ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (cell->at_boundary (f))
        {

          // we consider the center of the face and a vertex
          const Point<dim> cell_center = cell->center ();
          const Point<dim> face_center = cell->face (f)->center ();
          const Point<dim> face_vertex = cell->face (f)->vertex (0);
          // generate a vector parallel to the face
          const Tensor<1, dim> face_parallel = face_center - face_vertex;
          const Tensor<1, dim> face_normal = face_center - cell_center;
          unsigned int normal_axis = -1;
          for (unsigned int d = 0; d < dim; ++d)
            if (std::abs (face_parallel[d]) < eps)
              normal_axis = d;

          const double normal_component = face_normal[normal_axis];

          Assert(normal_axis == 0 or normal_axis == 1 or normal_axis == 2,
              ExcMessage("Wrong axis position."));
          Assert(std::abs (normal_component) >= 1.e-10,
              ExcMessage("Wrong direction. Normal can not be zero."))
          if (std::abs (normal_component) <= 1.e-10)
            std::cout << "Wrong direction. Normal can not be zero."
                      << std::endl;

          // How we do number this?
          // x- x+ y- y+ z- z+
          // 0  1  2  3  4  5
          unsigned int face_index = -1;
          if (normal_component <= 0.0)
            face_index = 2 * normal_axis;
          else
            face_index = 2 * normal_axis + 1;
          cell->face (f)->set_boundary_id (face_index);
        }
      }
    }
  }

  template <int dim>
  void
  ProbGeom<dim>::set_map_pin_index ()
  {
    // number of pins per row (x-axis)
    std::vector<unsigned int> pin_per_row;
    for (unsigned int i = 0; i < m_core.m_n_nodes[0]; ++i)
    {
      unsigned int max_i_size = 0;
      for (unsigned int j = 0; j < m_core.m_n_nodes[1]; ++j)
      {
        if (m_core.m_indices[0][j][i] != m_core.m_empty_assembly)
        {
          const unsigned int core_components =
              m_core.m_components[m_core.m_indices[0][j][i]];
          const std::vector<unsigned int> size_pins =
              mp_data.mp_geom.m_lattices.at(core_components).m_n_nodes;
          max_i_size = std::max (max_i_size, size_pins[0]);
        }
      }
      pin_per_row.push_back (max_i_size);
    }

    std::vector<unsigned int> cum_pin_per_row (pin_per_row.size ());
    for (unsigned int i = 1; i < pin_per_row.size (); ++i)
    {
      cum_pin_per_row[i] = cum_pin_per_row[i - 1] + pin_per_row[i - 1];
    }

    // pins per column (y-axis)
    std::vector<unsigned int> pin_per_col;
    for (unsigned int j = 0; j < m_core.m_n_nodes[1]; ++j)
    {
      unsigned int max_j_size = 0;
      for (unsigned int i = 0; i < m_core.m_n_nodes[0]; ++i)
      {
        if (m_core.m_indices[0][j][i] != m_core.m_empty_assembly)
        {
          const unsigned int core_components =
              m_core.m_components[m_core.m_indices[0][j][i]];
          const std::vector<unsigned int> size_pins =
              mp_data.mp_geom.m_lattices.at(core_components).m_n_nodes;
          max_j_size = std::max (max_j_size, size_pins[1]);
        }
      }
      pin_per_col.push_back (max_j_size);
    }

    std::vector<unsigned int> cum_pin_per_col (pin_per_col.size ());
    for (unsigned int i = 1; i < pin_per_col.size (); ++i)
    {
      cum_pin_per_col[i] = cum_pin_per_col[i - 1] + pin_per_col[i - 1];
    }

    // now we resize the map_pin_index, initializing it with a wrong value.
    const unsigned int n_cols = cum_pin_per_col.back () + pin_per_col.back ();
    const unsigned int n_rows = cum_pin_per_row.back () + pin_per_row.back ();

    map_pin_index.resize (n_cols,
        std::vector<unsigned int> (n_rows, m_core.m_empty_assembly));

    // now we assign the "user indices" associated with the particular pin
    unsigned int counter = 0;
    for (unsigned int j = 0; j < m_core.m_n_nodes[1]; ++j)
    {
      for (unsigned int i = 0; i < m_core.m_n_nodes[0]; ++i)
      {
        if (m_core.m_indices[0][j][i] != m_core.m_empty_assembly)
        {
          const unsigned int core_components =
              m_core.m_components[m_core.m_indices[0][j][i]];
          const std::vector<std::vector<unsigned int> > map_pins =
              mp_data.mp_geom.m_lattices.at(core_components).m_components;
          const std::vector<unsigned int> size_pins =
              mp_data.mp_geom.m_lattices.at(core_components).m_n_nodes;

          const unsigned int i_cum = cum_pin_per_row[i];
          const unsigned int j_cum = cum_pin_per_col[j];

          // std::cout << "filling non-empty assembly" << std::endl;
          for (unsigned int jj = 0; jj < size_pins[1]; ++jj)
          {
            for (unsigned int ii = 0; ii < size_pins[0]; ++ii)
            {
              if (map_pins[jj][ii] != m_core.m_empty_assembly)
              {
                map_pin_index[j_cum + jj][i_cum + ii] = counter;
                counter++;
              }
            }
          }
        }
      }
    }
  }

  template <int dim>
  double
  ProbGeom<dim>::area_correction_factor (const unsigned int n_ref)
  {
    const double pi = 3.14159265358979323846;
    const double pi_over_2_pow = pi / std::pow (2, n_ref + 1);
    const double correction_factor = std::sqrt (
        pi_over_2_pow / std::sin (pi_over_2_pow));
    if (dim == 1)
    {
      return 1.0;
    }
    else
    {
      return correction_factor;
    }
  }

  template <int dim>
  void
  ProbGeom<dim>::set_n_assemblies ()
  {
    /// @todo we count how many pins, but in the future we need something better than this loop.
    unsigned int tot_assemblies = 0;
    for (unsigned int j = 0; j < m_core.m_indices[0].size (); ++j)
    {
      for (unsigned int i = 0; i < m_core.m_indices[0][j].size (); ++i)
      {
        if (m_core.m_indices[0][j][i] != m_core.m_empty_assembly)
        {
          tot_assemblies += 1;
        }
      }
    }
    n_assemblies = tot_assemblies;
  }

  template <int dim>
  void
  ProbGeom<dim>::set_pins_per_assembly ()
  {
    // we count how many pins
    /// @todo we need something better than this loop?.
    const std::vector<unsigned int> n_latt = mp_data.mp_geom.m_core.m_n_nodes;
    for (unsigned int j = 0; j < n_latt[1]; ++j)
    {
      for (unsigned int i = 0; i < n_latt[0]; ++i)
      {
        if (m_core.m_indices[0][j][i] != m_core.m_empty_assembly)
        {
          unsigned int this_lattice = m_core.m_components[m_core.m_indices[0][j][i]];
          const std::vector<std::vector<unsigned int> > & map_pins =
              mp_data.mp_geom.m_lattices.at(this_lattice).m_components;
          const std::vector<unsigned int> & size_pins =
              mp_data.mp_geom.m_lattices.at(this_lattice).m_n_nodes;
          unsigned int tot_pins = 0;
          for (unsigned int ii = 0; ii < size_pins[1]; ++ii)
          {
            for (unsigned int jj = 0; jj < size_pins[0]; ++jj)
            {
              if (map_pins[ii][jj] != m_core.m_empty_assembly)
              {
                ++tot_pins;
              }
            }
          }
          pins_per_assembly.push_back (tot_pins);
        }
      }
    }
    Assert(pins_per_assembly.size () == n_assemblies,
        ExcMessage("Wrong number of assemplies counted."));
  }

  template <int dim>
  void
  ProbGeom<dim>::set_n_pins ()
  {
    unsigned int tot_pins = 0;
    for (unsigned int i = 0; i < n_assemblies; ++i)
    {
      tot_pins += pins_per_assembly[i];
    }
    n_pins = tot_pins;
  }

  /**
   * @brief printing the mesh the used for the calculation
   * @param outputeps
   * @todo Add an specialization to print the 1D geometry as a 2D slab
   */
  template <int dim>
  void
  ProbGeom<dim>::print_mesh_eps (const std::string & outputeps) const
  {
    // we do not print the eps for 3d because it is too messy
    //if (dim == 2)
    {
      std::ofstream eps_out (outputeps);
      GridOut gridout;
      gridout.write_eps (m_tria, eps_out);
    }

  }

  /**
   * @brief printing the mesh the used for the calculation
   * @param outputvtk
   * @todo Add an specialization to print the 1D geometry as a 2D slab
   */
  template <int dim>
  void
  ProbGeom<dim>::print_mesh_vtk (const std::string & outputvtk) const
  {
    // vtk
    {
      std::ofstream vtk_out (outputvtk);
      DoFHandler<dim, dim> dof_handler (m_tria);
      static const FE_DGQ<dim> fe (0);
      dof_handler.distribute_dofs (fe);
      DataOut<dim, DoFHandler<dim, dim> > data_out;
      data_out.attach_dof_handler (dof_handler);
      Vector<double> mat_id_cells (m_tria.n_active_cells ());
      Vector<double> user_index_cells (m_tria.n_active_cells ());
      unsigned int j = 0;
      for (typename Triangulation<dim, dim>::active_cell_iterator cell =
          m_tria.begin_active (); cell != m_tria.end (); ++cell, ++j)
      {
        const unsigned int mat_id = cell->material_id ();
        mat_id_cells[j] = static_cast<double> (mat_id);
        user_index_cells[j] = static_cast<double> (cell->user_index ());
      }
      data_out.add_data_vector (mat_id_cells, "material",
          DataOut<dim>::type_cell_data);
      data_out.add_data_vector (user_index_cells, "UserIndex",
          DataOut<dim>::type_cell_data);
      data_out.build_patches ();
      data_out.write_vtk (vtk_out);
    }
  }

  template <>
  void
  ProbGeom<1>::print_mesh_eps (const std::string & outputeps) const
  {
    // We calculate how long is the reactor
    double total_length = 0, a = 0, b = 0;
    for (typename Triangulation<1, 1>::active_cell_iterator cell =
        m_tria.begin_active (); cell != m_tria.end (); ++cell)
    {
      for (unsigned int v = 0; v < GeometryInfo<1>::vertices_per_cell; ++v)
      {
        if (cell->vertex (v)[0] < a)
          a = cell->vertex (v)[0];
        if (cell->vertex (v)[0] > b)
          b = cell->vertex (v)[0];
      }
    }
    total_length = b - a;

    // We flatten the triangulation in order to extrude it later
    Triangulation<1, 1> tria_tmp;
    GridGenerator::flatten_triangulation (m_tria, tria_tmp);
    // Now we extrude the triangulation. For this, the basic settings is
    // 2 slices and the height using an standard proportion
    const unsigned int n_slices = 2;
    const double height = total_length / ((1. + std::sqrt (5)) / 2.);
    Triangulation<2, 2> tria2d;
    GridGenerator::extrude_triangulation (tria_tmp, n_slices, height, tria2d);

    // we can print the eps file
    {
      std::ofstream eps_out (outputeps);
      GridOut gridout;
      gridout.write_eps (tria2d, eps_out);
    }



  }

  template <>
  void
  ProbGeom<1>::print_mesh_vtk (const std::string & outputvtk) const
  {
    // We calculate how long is the reactor
    double total_length = 0, a = 0, b = 0;
    for (typename Triangulation<1, 1>::active_cell_iterator cell =
        m_tria.begin_active (); cell != m_tria.end (); ++cell)
    {
      for (unsigned int v = 0; v < GeometryInfo<1>::vertices_per_cell; ++v)
      {
        if (cell->vertex (v)[0] < a)
          a = cell->vertex (v)[0];
        if (cell->vertex (v)[0] > b)
          b = cell->vertex (v)[0];
      }
    }
    total_length = b - a;

    // We flatten the triangulation in order to extrude it later
    Triangulation<1, 1> tria_tmp;
    GridGenerator::flatten_triangulation (m_tria, tria_tmp);
    // Now we extrude the triangulation. For this, the basic settings is
    // 2 slices and the height using an standard proportion
    const unsigned int n_slices = 2;
    const double height = total_length / ((1. + std::sqrt (5)) / 2.);
    Triangulation<2, 2> tria2d;
    GridGenerator::extrude_triangulation (tria_tmp, n_slices, height, tria2d);

    // Now we print the vtk file containing the material id and some recreated
    // user indices.
    {
      std::ofstream vtk_out (outputvtk);
      DoFHandler<2, 2> dof_handler (tria2d);
      static const FE_DGQ<2> fe (0);
      dof_handler.distribute_dofs (fe);
      DataOut<2, DoFHandler<2, 2> > data_out;
      data_out.attach_dof_handler (dof_handler);
      Vector<double> mat_id_cells (tria2d.n_active_cells ());
      Vector<double> user_index_cells (tria2d.n_active_cells ());
      unsigned int j = 0;
      for (typename Triangulation<2, 2>::active_cell_iterator cell =
          tria2d.begin_active (); cell != tria2d.end (); ++cell, ++j)
      {
        const unsigned int mat_id = cell->material_id ();
        mat_id_cells[j] = static_cast<double> (mat_id);
        cell->set_user_index (j);
        user_index_cells[j] = static_cast<double> (cell->user_index ());
      }
      data_out.add_data_vector (mat_id_cells, "material",
          DataOut<2>::type_cell_data);
      data_out.add_data_vector (user_index_cells, "UserIndex",
          DataOut<2>::type_cell_data);
      data_out.build_patches ();
      data_out.write_vtk (vtk_out);
    }

  }

  template <>
  void
  ProbGeom<1>::set_manifold (Triangulation<1, 1> & /* tria */,
                             std::vector<RealPin<1> > & /* real_pins */,
                             const unsigned int /* first_manifold */)
  {
    Assert(false, ExcMessage("Not implemented."));
  }

  template <>
  void
  ProbGeom<2>::set_manifold (Triangulation<2, 2> & tria,
                             std::vector<RealPin<2> > & real_pins,
                             const unsigned int first_manifold)
  {

    // set the different identifiers to the mesh, in order to attach
    // different manifolds to different identifiers
    assign_pin_manifolds (tria, real_pins, first_manifold);

    // we associate the internal manifolds with boundary entities
    //std::vector<SphericalManifold<2, 2> > manifold;
    for (unsigned int i = 0; i < real_pins.size (); ++i)
    {
      //SphericalManifold<2, 2> manifold_aux (real_pins[i].get_center ());
      //manifold.push_back (manifold_aux);
      std::shared_ptr<Manifold<2> > manifold_aux (
          new SphericalManifold<2, 2> (real_pins[i].get_center ()));
      m_manifold.push_back (manifold_aux);
    }
    for (unsigned int i = 0; i < real_pins.size (); ++i)
    {
      if (real_pins[i].get_type () == PinType::PIN or real_pins[i].get_type ()
          == PinType::PINNEW)
      {
        unsigned int manifold_index = first_manifold + i;
        //tria.set_manifold (manifold_index, manifold[i]);
        tria.set_manifold (manifold_index, *m_manifold[i]);
      }
    }
  }

  template <>
  void
  ProbGeom<3>::set_manifold (Triangulation<3, 3> & tria,
                             std::vector<RealPin<3> > & real_pins,
                             const unsigned int first_manifold)
  {
    // set the different identifiers to the mesh, in order to attach
    // different manifolds to different identifiers
    assign_pin_manifolds (tria, real_pins, first_manifold);

    // we associate the internal manifolds with boundary entities
    //std::vector<CylindricalManifold<3> > manifold;
    for (unsigned int i = 0; i < real_pins.size (); ++i)
    {
      const Point<3> tmp_center (real_pins[i].get_center ());
      const Point<3> tmp_dir (0, 0, 1);
      //CylindricalManifold<3> manifold_aux (tmp_dir, tmp_center);
      //manifold.push_back (manifold_aux);
      std::shared_ptr<Manifold<3> > manifold_aux (
          new CylindricalManifold<3> (tmp_dir, tmp_center));
      m_manifold.push_back (manifold_aux);
    }
    for (unsigned int i = 0; i < real_pins.size (); ++i)
    {
      if (real_pins[i].get_type () == PinType::PIN or real_pins[i].get_type ()
          == PinType::PINNEW)
      {
        unsigned int manifold_index = first_manifold + i;
        //tria.set_manifold (manifold_index, manifold[i]);
        tria.set_manifold (manifold_index, *m_manifold[i]);
      }
    }
  }

  /**
   * Constructor
   */
  template <int dim>
  ProbGeom<dim>::ProbGeom (const Input & data)
      : mp_data (data),
        m_core (mp_data.mp_geom.m_core),
        m_mesh_regularity(0)
  {
    // do we talk?
    const bool verbose = false;

    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart (); // We restart the timer

    // message
    if (verbose)
      log_print_text ("Generating geometry", 1);

    // get the number of refinements
    const unsigned int n_ref = mp_data.mp_settings.get_n_ref ();

    // calculate the correction factor for the pin
    const double correction_factor = area_correction_factor (n_ref);

    // Real core generation
    RealCore<dim> real_core (mp_data, correction_factor);

    // Generate the "reactor" information for the geometry
    n_assemblies = real_core.get_n_assemblies ();
    n_pins = real_core.get_n_pins ();
    pin2assembly = real_core.get_pins_to_assembly ();
    pins_per_assembly = real_core.get_pins_per_assembly ();

    // We build the mesh
    real_core.build_mesh (m_tria);

    // attach manifolds
    if (dim > 1)
    {
      const unsigned int first_manifold = 1024;
      set_manifold (m_tria, real_core.m_real_pins, first_manifold);
    }

    // refine the mesh
    refine_assembly (m_tria, real_core.m_real_pins, n_ref);

    // Information linking the cells of the triangulation with reactor info.
    /// @todo Here we should assign the indices before refinement, and then use
    /// recursively_set_user_index() at one level up (should be cheaper).
    cell2pin = assign_cell_index (m_tria, real_core.m_real_pins);

    // set the boundary indicators
    set_boundary_indicators ();
    set_map_pin_index ();

    // Print the mesh if required
    if (mp_data.mp_settings.get_print_mesh_flag ())
    {
      print_mesh_eps (mp_data.mp_settings.get_file_mesheps ());
      print_mesh_vtk (mp_data.mp_settings.get_file_meshvtk ());
    }

    /// We now check the regularity of the mesh (1 for all cells equal, 0 otherwise)
    this->set_mesh_regularity();
    //deallog << "the mesh regularity is " << this->get_mesh_regularity() << std::endl;

    // Stop the timer and report the time
    timer.stop ();
    StatsTimer::instance ().add ("Generating geometry", timer.wall_time ());
    // Reprot the memory
    StatsMemory::instance ().add ("Geometry", memory_consumption ());
  }

  /// For the moment we use all cells the same == 1, otherwise == 0
  template <int dim>
  void
  ProbGeom<dim>::set_mesh_regularity ()
  {
    m_mesh_regularity = 1;
    const typename Triangulation<dim, dim>::active_cell_iterator first_cell = m_tria.begin_active ();
    for (typename Triangulation<dim, dim>::active_cell_iterator cell = m_tria.begin_active (); cell != m_tria.end (); ++cell)
    {
      if (cell == m_tria.begin_active())
      {
        continue;
      }
      // just one different element means no regularity
      if (not cell->is_translation_of(first_cell))
      {
        m_mesh_regularity = 0;
        break;
      }
    }
  }

  template <int dim>
  double
  ProbGeom<dim>::memory_consumption () const
  {
    double memory_consumption = 0;
    memory_consumption += m_tria.memory_consumption ();
    memory_consumption += sizeof(map_pin_index);
    return memory_consumption;
  }

  template class ProbGeom<1> ;
  template class ProbGeom<2> ;
  template class ProbGeom<3> ;

} // end of name space Forest

