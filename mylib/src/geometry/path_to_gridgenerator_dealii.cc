/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   geometry/path_to_gridgenerator_dealii.cc
 * @brief  Implementation of grid generator functions
 */

#include "geometry/path_to_gridgenerator_dealii.h"

#include "deal.II/base/exceptions.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/point.h"
#include "deal.II/distributed/tria.h"
#include "deal.II/grid/grid_reordering.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_generator.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace dealii
{
  namespace GridGenerator
  {

    void
    extrude_triangulation (const Triangulation<1, 1> &input,
                           const unsigned int n_slices,
                           const double height,
                           Triangulation<2, 2> &result)
    {
      Assert(input.n_levels() == 1,
          ExcMessage ("The input triangulation must be a coarse mesh, i.e., it must " "not have been refined."));
      Assert(result.n_cells()==0,
          ExcMessage("The output triangulation object needs to be empty."));
      Assert(height>0,
          ExcMessage("The given height for extrusion must be positive."));
      Assert(n_slices>=1,
          ExcMessage("The number of slices for extrusion must be at least 2."));

      std::vector<Point<2> > points (n_slices * input.n_vertices ());
      std::vector<CellData<2> > cells;
      cells.reserve ((n_slices - 1) * input.n_active_cells ());

      // copy the array of points as many times as there will be slices,
      // one slice at a time
      for (unsigned int slice = 0; slice < n_slices; ++slice)
      {
        for (unsigned int i = 0; i < input.n_vertices (); ++i)
        {
          const Point<1> &v = input.get_vertices ()[i];
          points[slice * input.n_vertices () + i] (0) = v (0);
          points[slice * input.n_vertices () + i] (1) = height * slice
              / (n_slices - 1);
        }
      }

      // then create the cells of each of the slices, one stack at a
      // time
      for (Triangulation<1, 1>::cell_iterator cell = input.begin ();
          cell != input.end (); ++cell)
      {
        for (unsigned int slice = 0; slice < n_slices - 1; ++slice)
        {
          CellData<2> this_cell;
          for (unsigned int v = 0; v < GeometryInfo<1>::vertices_per_cell; ++v)
          {
            this_cell.vertices[v] = cell->vertex_index (v)
                + slice * input.n_vertices ();
            this_cell.vertices[v + GeometryInfo<1>::vertices_per_cell] =
                cell->vertex_index (v) + (slice + 1) * input.n_vertices ();
          }

          this_cell.material_id = cell->material_id ();
          cells.push_back (this_cell);
        }
      }

      // use all of this to finally create the extruded 2d triangulation
      result.create_triangulation (points, cells, SubCellData ());
    }

    /**
     * @details Given an input triangulation in_tria, this function makes
     * a new flat triangulation out_tria which contains a single
     * level with all active cells of the input triangulation.
     *
     * All information about cell manifold_ids and material ids
     * are copied from one triangulation to the other, and only
     * the boundary manifold_ids and boundary_ids are copied over
     * from the faces of @p in_tria to the faces of @p out_tria. If
     * you need to specify manifold ids on interior faces, they have
     * to be specified manually after the triangulation is created.
     */
    template <int dim, int spacedim1, int spacedim2>
    void
    flatten_triangulation (const Triangulation<dim, spacedim1> &in_tria,
                           Triangulation<dim, spacedim2> &out_tria)
    {
      // we don't use it so we comment it, thus skipping the compiling warning
      const parallel::distributed::Triangulation<dim, spacedim1> *pt =
          dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim1> *> (&in_tria);

      Assert(pt == NULL,
          ExcMessage("Cannot use this function on parallel::distributed::Triangulation."));
      assert(pt == NULL);

      std::vector<Point<spacedim2> > v;
      std::vector<CellData<dim> > cells;
      SubCellData subcelldata;

      const unsigned int spacedim = std::min (spacedim1, spacedim2);
      const std::vector<Point<spacedim1> > &in_vertices =
          in_tria.get_vertices ();

      v.resize (in_vertices.size ());
      for (unsigned int i = 0; i < in_vertices.size (); ++i)
        for (unsigned int d = 0; d < spacedim; ++d)
          v[i][d] = in_vertices[i][d];

      cells.resize (in_tria.n_active_cells ());
      typename Triangulation<dim, spacedim1>::active_cell_iterator cell =
          in_tria.begin_active (), endc = in_tria.end ();

      for (unsigned int id = 0; cell != endc; ++cell, ++id)
      {
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
          cells[id].vertices[i] = cell->vertex_index (i);
        cells[id].material_id = cell->material_id ();
        cells[id].manifold_id = cell->manifold_id ();
      }

      if (dim > 1)
      {
        typename Triangulation<dim, spacedim1>::active_face_iterator face =
            in_tria.begin_active_face (), endf = in_tria.end_face ();

        // Face counter for both dim == 2 and dim == 3
        unsigned int f = 0;
        switch (dim)
        {
          case 2:
          {
            subcelldata.boundary_lines.resize (in_tria.n_active_faces ());
            for (; face != endf; ++face)
              if (face->at_boundary ())
              {
                for (unsigned int i = 0;
                    i < GeometryInfo<dim>::vertices_per_face; ++i)
                  subcelldata.boundary_lines[f].vertices[i] =
                      face->vertex_index (i);
                subcelldata.boundary_lines[f].boundary_id =
                    face->boundary_id ();
                subcelldata.boundary_lines[f].manifold_id =
                    face->manifold_id ();
                ++f;
              }
            subcelldata.boundary_lines.resize (f);
          }
            break;
          case 3:
          {
            subcelldata.boundary_quads.resize (in_tria.n_active_faces ());
            for (; face != endf; ++face)
              if (face->at_boundary ())
              {
                for (unsigned int i = 0;
                    i < GeometryInfo<dim>::vertices_per_face; ++i)
                  subcelldata.boundary_quads[f].vertices[i] =
                      face->vertex_index (i);
                subcelldata.boundary_quads[f].boundary_id =
                    face->boundary_id ();
                subcelldata.boundary_quads[f].manifold_id =
                    face->manifold_id ();
                ++f;
              }
            subcelldata.boundary_quads.resize (f);
          }
            break;
          default:
          {
            Assert(false, ExcInternalError());
            break;
          }
        }
      }
      out_tria.create_triangulation (v, cells, subcelldata);
    }

    /**
     * Overload of the function merge_triangulation to deal with
     * a vector of triangulations (not just two).
     */
    template <int dim, int spacedim>
    void
    merge_triangulations (const std::vector<Triangulation<dim, spacedim> > &tria,
                          Triangulation<dim, spacedim> &result)
    {

      const unsigned int n_tri = tria.size ();

      // Check that the number of triangulation is at least one
      Assert(n_tri >= 1, ExcMessage ("n_tri should be at least one"));

      if (n_tri == 1)
      {
        result.clear ();
        result.copy_triangulation (tria[0]);
        return;
      }

      // Check that the triangulations are coarse meshes
      for (unsigned int i = 0; i < n_tri; ++i)
      {
        Assert(tria[i].n_levels() == 1,
            ExcMessage ("The input triangulations must be coarse meshes."));
      }

      // get the union of the set of vertices
      std::vector<Point<spacedim> > vertices = tria[0].get_vertices ();
      for (unsigned int i = 1; i < n_tri; ++i)
      {
        vertices.insert (vertices.end (), tria[i].get_vertices ().begin (),
            tria[i].get_vertices ().end ());
      }

      // How many cells
      unsigned int n_cells = 0;
      for (unsigned int i = 0; i < n_tri; ++i)
        n_cells += tria[i].n_cells ();

      std::vector<CellData<dim> > cells;
      cells.reserve (n_cells);

      // now form the union of the set of cells. note that we have to
      // translate the vertex indices
      unsigned int n_vertices = 0;
      for (unsigned int i = 0; i < n_tri; ++i)
      {
        for (typename Triangulation<dim, spacedim>::cell_iterator cell =
            tria[i].begin (); cell != tria[i].end (); ++cell)
        {
          CellData<dim> this_cell;
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
              ++v)
            this_cell.vertices[v] = cell->vertex_index (v) + n_vertices;
          this_cell.material_id = cell->material_id ();
          cells.push_back (this_cell);
        }
        n_vertices += tria[i].n_vertices ();
      }

      // we do the loop because we have more than two triangulations, and
      // delete_duplicated_vertices only works for 2 triangulations.
      // for (unsigned int i = 0; i < n_tri-1; ++i)
      // if we assume that only 2^dim cells can share a vertex, using
      // delete_duplicated_vertices just "dim" times should be enough
      for (unsigned int i = 0; i < dim; ++i)
      {
        // throw out duplicated vertices from the two meshes, reorder vertices as
        // necessary and create the triangulation
        SubCellData subcell_data;
        std::vector<unsigned int> considered_vertices;
        GridTools::delete_duplicated_vertices (vertices, cells, subcell_data,
            considered_vertices);
      }
      // One more time
      SubCellData subcell_data;
      std::vector<unsigned int> considered_vertices;
      GridTools::delete_duplicated_vertices (vertices, cells, subcell_data,
          considered_vertices);

      //for (unsigned int i = 0; i<vertices.size(); ++i)
      //  std::cout << "vertex["<<i<<"] = " << vertices[i] << std::endl;

      // reorder the cells to ensure that they satisfy the convention for
      // edge and face directions
      GridReordering<dim, spacedim>::reorder_cells (cells, true);
      result.clear ();
      result.create_triangulation (vertices, cells, subcell_data);
    }

    template <int dim>
    void
    pin_cell (Triangulation<dim, dim> &,
              const Point<dim> &,
              const Point<dim> &,
              const double,
              const std::vector<unsigned int> &)
    {
      Assert(false, ExcMessage ("Not implemented."));
    }

    template <>
    void
    pin_cell (Triangulation<2, 2> & tria,
              const Point<2> & center,
              const Point<2> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id)
    {
      //const unsigned int n_verts = 8;
      const unsigned int n_cells = 5;

      const double px = pitch[0] / 2.0;
      const double py = pitch[1] / 2.0;
      const double r = pin_radius * std::sqrt (2) / 2;

      const double xx = center (0);
      const double yy = center (1);

      /**
       *  Defining special vertices manually
       *
       *  2-----------3
       *  | \       / |
       *  |  6-----7  |
       *  |  |     |  |
       *  |  |     |  |
       *  |  4-----5  |
       *  | /       \ |
       *  0-----------1
       */
      const std::vector<Point<2> > vertices =
      { Point<2> (xx - px, yy - py), /**<  v0 box */
        Point<2> (xx + px, yy - py), /**<  v1 box */
        Point<2> (xx - px, yy + py), /**<  v2 box */
        Point<2> (xx + px, yy + py), /**<  v3 box */
        Point<2> (xx - r, yy - r), /**<  v4 circle */
        Point<2> (xx + r, yy - r), /**<  v5 circle */
        Point<2> (xx - r, yy + r), /**<  v6 circle */
        Point<2> (xx + r, yy + r) /**<  v7 circle */
      };

      /**
       *  Define connections following the deal.II standards
       *       3
       *    2-->--3
       *    |     |
       *   0^     ^1
       *    |     |
       *    0-->--1
       *       2
       *
       * This cell is {0, 1, 2, 3}
       */
      const int cell_vertices[n_cells][GeometryInfo<2>::vertices_per_cell] =
      {
      { 0, 1, 4, 5 },
        { 0, 4, 2, 6 },
        { 4, 5, 6, 7 },
        { 1, 3, 5, 7 },
        { 2, 6, 3, 7 } };

      std::vector<CellData<2> > cells (n_cells, CellData<2> ());

      // creating arrays with vertex info
      for (unsigned int i = 0; i < n_cells; ++i)
        for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];

      for (unsigned int i = 0; i < n_cells; ++i)
      {
        // i == 2 is the cell in the center (fuel), and the others are water
        if (i == 2)
          cells[i].material_id = mat_id[1];
        else
          cells[i].material_id = mat_id[0];
      }

      // Create the triangulation with no boundary information
      tria.create_triangulation (vertices, cells, SubCellData ());
    }

    template <>
    void
    pin_cell (Triangulation<1, 1> & tria,
              const Point<1> & center,
              const Point<1> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id)
    {
      //const unsigned int n_verts = 4;
      const unsigned int n_cells = 3;

      const double px = pitch[0] / 2.0;
      const double r = pin_radius;

      const double xx = center (0);

      /**
       *  Defining special vertices manually
       *
       *  0---2---3---1
       */
      const std::vector<Point<1> > vertices =
      { Point<1> (xx - px), /**<  v0 box */
        Point<1> (xx + px), /**<  v1 box */
        Point<1> (xx - r),/**<  v2 circle */
        Point<1> (xx + r) /**<  v3 circle */
      };

      /**
       *  Define connections following the deal.II standards
       *       1
       *    0-->--1
       *
       * This cell is {0, 1, 2, 3}
       */
      const int cell_vertices[n_cells][GeometryInfo<1>::vertices_per_cell] =
      {
      { 0, 2 },
        { 2, 3 },
        { 3, 1 } };

      std::vector<CellData<1> > cells (n_cells, CellData<1> ());

      // creating arrays with vertex info
      for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < GeometryInfo<1>::vertices_per_cell; ++j)
        {
          cells[i].vertices[j] = cell_vertices[i][j];
        }

        // i == 1 is the cell in the center (fuel), and the others are water
        if (i == 1)
          cells[i].material_id = mat_id[1];
        else
          cells[i].material_id = mat_id[0];
      }

      // Create the triangulation with no boundary information
      tria.create_triangulation (vertices, cells, SubCellData ());
    }

    template <>
    void
    pin_cell (Triangulation<3, 3> & tria,
              const Point<3> & center,
              const Point<3> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id)
    {

      const Point<2> center_2d (center (0), center (1));
      const Point<2> pitch_2d (pitch (0), pitch (1));

      Triangulation<2, 2> tria_2d;
      pin_cell (tria_2d, center_2d, pitch_2d, pin_radius, mat_id);

      unsigned int n_nodes_z = 1;
      double step_z = pitch(2);
      Assert(n_nodes_z > 0, ExcMessage("n_nodes_z must be minimum 1"));

      //std::cout << "Now we extrude the mesh to have a 3d assembly\n";
      GridGenerator::extrude_triangulation (tria_2d, n_nodes_z + 1, step_z,
          tria);
    }

    template <int dim>
    void
    pin_cell2 (Triangulation<dim, dim> &,
               const Point<dim> &,
               const Point<dim> &,
               const double,
               const std::vector<unsigned int> &)
    {
      Assert(false, ExcMessage ("Not implemented."));
    }

    template <>
    void
    pin_cell2 (Triangulation<1, 1> & tria,
               const Point<1> & center,
               const Point<1> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id)
    {
      pin_cell<1> (tria, center, pitch, pin_radius, mat_id);
    }

    template <>
    void
    pin_cell2 (Triangulation<2, 2> & tria,
               const Point<2> & center,
               const Point<2> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id)
    {
      // const unsigned int n_verts = 12;
      const unsigned int n_cells = 9;

      const double px = pitch[0] / 2.0;
      const double py = pitch[1] / 2.0;
      const double r = pin_radius * std::sqrt (2) / 2;
      const double r2 = r / 3.0;

      const double xx = center (0);
      const double yy = center (1);

      /**
       *  Defining special vertices manually
       *
       *  2-------------3
       *  | \         / |
       *  |  6-------7  |
       *  |  | \   / |  |
       *  |  | 10-11 |  |
       *  |  |  8-9  |  |
       *  |  | /   \ |  |
       *  |  4-------5  |
       *  | /         \ |
       *  0-------------1
       */
      const std::vector<Point<2> > vertices =
      { Point<2> (xx - px, yy - py), /**<  v0 box */
        Point<2> (xx + px, yy - py), /**<  v0 box */
        Point<2> (xx - px, yy + py), /**<  v0 box */
        Point<2> (xx + px, yy + py), /**<  v0 box */
        // circle
        Point<2> (xx - r, yy - r), /**<  v0 circle */
        Point<2> (xx + r, yy - r), /**<  v0 circle */
        Point<2> (xx - r, yy + r), /**<  v0 circle */
        Point<2> (xx + r, yy + r), /**<  v0 circle */
        // circle2
        Point<2> (xx - r2, yy - r2), /**<  v0 circle */
        Point<2> (xx + r2, yy - r2), /**<  v0 circle */
        Point<2> (xx - r2, yy + r2), /**<  v0 circle */
        Point<2> (xx + r2, yy + r2) /**<  v0 circle */};

      /**
       *  Define connections following the deal.II standards
       *       3
       *    2-->--3
       *    |     |
       *   0^     ^1
       *    |     |
       *    0-->--1
       *       2
       *
       * This cell is {0, 1, 2, 3}
       */
      const int cell_vertices[n_cells][GeometryInfo<2>::vertices_per_cell] =
      {
      { 0, 1, 4, 5 },
        { 0, 4, 2, 6 },
        { 2, 6, 3, 7 },
        { 1, 3, 5, 7 },
        { 4, 5, 8, 9 },
        { 4, 8, 6, 10 },
        { 6, 10, 7, 11 },
        { 5, 7, 9, 11 },
        { 8, 9, 10, 11 } };

      std::vector<CellData<2> > cells (n_cells, CellData<2> ());

      // creating arrays with vertex info
      for (unsigned int i = 0; i < n_cells; ++i)
        for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];

      for (unsigned int i = 0; i < n_cells; ++i)
      {
        // i > 3 is the cell in the center (fuel), and the others are water
        if (i > 3)
          cells[i].material_id = mat_id[1];
        else
          cells[i].material_id = mat_id[0];
      }

      // Create the triangulation with no boundary information
      tria.create_triangulation (vertices, cells, SubCellData ());
    }

    template <>
    void
    pin_cell2 (Triangulation<3, 3> & tria,
               const Point<3> & center,
               const Point<3> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id)
    {

      const Point<2> center_2d (center (0), center (1));
      const Point<2> pitch_2d (pitch (0), pitch (1));

      Triangulation<2, 2> tria_2d;
      pin_cell2 (tria_2d, center_2d, pitch_2d, pin_radius, mat_id);

      unsigned int n_nodes_z = 1;
      double step_z = pitch(2);
      Assert(n_nodes_z > 0, ExcMessage("n_nodes_z must be minimum 1"));

      //std::cout << "Now we extrude the mesh to have a 3d assembly\n";
      GridGenerator::extrude_triangulation (tria_2d, n_nodes_z + 1, step_z,
          tria);
    }

    template <int dim>
    void
    pin_box (Triangulation<dim> &,
             const Point<dim> &,
             const Point<dim> &,
             const std::vector<unsigned int> &)
    {
      Assert(false, ExcMessage ("Not implemented."));
    }

    template <>
    void
    pin_box (Triangulation<2, 2> & tria,
             const Point<2> & center,
             const Point<2> & pitch,
             const std::vector<unsigned int> & mat_id)
    {
      //const unsigned int n_verts = 4;
      const unsigned int n_cells = 1;

      const double px = pitch[0] / 2.0;
      const double py = pitch[1] / 2.0;

      const double xx = center (0);
      const double yy = center (1);

      /**
       *  Defining special vertices manually
       *
       *  2---3
       *  |   |
       *  0---1
       */

      const std::vector<Point<2> > vertices =
      { Point<2> (xx - px, yy - py), /**<  v0 box */
        Point<2> (xx + px, yy - py), /**<  v1 box */
        Point<2> (xx - px, yy + py), /**<  v2 box */
        Point<2> (xx + px, yy + py) /**<  v3 box */
      };

      /**
       *  Define connections following the deal.II standards
       *       3
       *    2-->--3
       *    |     |
       *   0^     ^1
       *    |     |
       *    0-->--1
       *       2
       *
       * This cell is {0, 1, 2, 3}
       */
      const int cell_vertices[n_cells][GeometryInfo<2>::vertices_per_cell] =
      {
      { 0, 1, 2, 3 } };

      std::vector<CellData<2> > cells (n_cells, CellData<2> ());

      // creating arrays with vertex info
      for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
        {
          cells[i].vertices[j] = cell_vertices[i][j];
        }
        cells[i].material_id = mat_id[0];
      }

      // Create the triangulation with no boundary information
      tria.create_triangulation (vertices, cells, SubCellData ());
    }

    template <>
    void
    pin_box (Triangulation<1, 1> & tria,
             const Point<1> & center,
             const Point<1> & pitch,
             const std::vector<unsigned int> & mat_id)
    {
      //const unsigned int n_verts = 2;
      const unsigned int n_cells = 1;

      const double px = pitch[0] / 2.0;

      const double xx = center (0);

      /**
       *  Defining special vertices manually
       *
       *  0---1
       */
      const std::vector<Point<1> > vertices =
      { Point<1> (xx - px), /**<  v0 box */
        Point<1> (xx + px) /**<  v1 box */
      };

      /**
       *  Define connections following the deal.II standards
       *       1
       *    0-->--1
       *
       * This cell is {0, 1, 2, 3}
       */
      const int cell_vertices[n_cells][GeometryInfo<1>::vertices_per_cell] =
      {
      { 0, 1 } };

      std::vector<CellData<1> > cells (n_cells, CellData<1> ());

      // creating arrays with vertex info
      for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < GeometryInfo<1>::vertices_per_cell; ++j)
        {
          cells[i].vertices[j] = cell_vertices[i][j];
        }
        cells[i].material_id = mat_id[0];
      }

      // Create the triangulation with no boundary information
      tria.create_triangulation (vertices, cells, SubCellData ());
    }

    template <>
    void
    pin_box (Triangulation<3, 3> & tria,
             const Point<3> & center,
             const Point<3> & pitch,
             const std::vector<unsigned int> & mat_id)
    {

      const Point<2> center_2d (center (0), center (1));
      const Point<2> pitch_2d (pitch (0), pitch (1));

      Triangulation<2, 2> tria_2d;
      pin_box (tria_2d, center_2d, pitch_2d, mat_id);

      unsigned int n_nodes_z = 1;
      double step_z = pitch(2);
      Assert(n_nodes_z > 0, ExcMessage("n_nodes_z must be minimum 1"));

      //std::cout << "Now we extrude the mesh to have a 3d assembly\n";
      GridGenerator::extrude_triangulation (tria_2d, n_nodes_z + 1, step_z,
          tria);
    }

    template
    void
    flatten_triangulation (const Triangulation<1, 1> &in_tria,
                           Triangulation<1, 1> &out_tria);

    template
    void
    flatten_triangulation (const Triangulation<2, 2> &in_tria,
                           Triangulation<2, 2> &out_tria);

    template
    void
    flatten_triangulation (const Triangulation<3, 3> &in_tria,
                           Triangulation<3, 3> &out_tria);

    template
    void
    merge_triangulations (const std::vector<Triangulation<1, 1> > &tria,
                          Triangulation<1, 1> &result);

    template
    void
    merge_triangulations (const std::vector<Triangulation<2, 2> > &tria,
                          Triangulation<2, 2> &result);

    template
    void
    merge_triangulations (const std::vector<Triangulation<3, 3> > &tria,
                          Triangulation<3, 3> &result);

    template
    void
    pin_cell (Triangulation<1, 1> & tria,
              const Point<1> & center,
              const Point<1> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id);

    template
    void
    pin_cell (Triangulation<2, 2> & tria,
              const Point<2> & center,
              const Point<2> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id);

    template
    void
    pin_cell (Triangulation<3, 3> & tria,
              const Point<3> & center,
              const Point<3> & pitch,
              const double pin_radius,
              const std::vector<unsigned int> & mat_id);

    template
    void
    pin_cell2 (Triangulation<1, 1> & tria,
               const Point<1> & center,
               const Point<1> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id);

    template
    void
    pin_cell2 (Triangulation<2, 2> & tria,
               const Point<2> & center,
               const Point<2> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id);

    template
    void
    pin_cell2 (Triangulation<3, 3> & tria,
               const Point<3> & center,
               const Point<3> & pitch,
               const double pin_radius,
               const std::vector<unsigned int> & mat_id);

    template
    void
    pin_box (Triangulation<1, 1> & tria,
             const Point<1> & center,
             const Point<1> & pitch,
             const std::vector<unsigned int> & mat_id);

    template
    void
    pin_box (Triangulation<2, 2> & tria,
             const Point<2> & center,
             const Point<2> & pitch,
             const std::vector<unsigned int> & mat_id);

    template
    void
    pin_box (Triangulation<3, 3> & tria,
             const Point<3> & center,
             const Point<3> & pitch,
             const std::vector<unsigned int> & mat_id);

  }
}

