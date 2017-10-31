/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/downstream_new.cc
 * @brief  Implementation of some renumbering functions
 */

#include "neutronics/downstream_new.h"

#include "deal.II/base/config.h"        // for DEAL_II_VERSION_GTE
#if DEAL_II_VERSION_GTE(8,4,0)
#include "deal.II/base/tensor.h"        // for Tensor
#endif
#include "deal.II/base/point.h"         // for Point
#include "deal.II/base/geometry_info.h"  // for GeometryInfo
#include "deal.II/base/quadrature_lib.h"  // for QGauss
#include "deal.II/base/types.h"         // for global_dof_index
#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/dofs/dof_renumbering.h"  // for compute_cell_wise
#include "deal.II/fe/fe_values.h"
#include "deal.II/fe/fe_update_flags.h"  // for operator|, UpdateFlags, etc

#include <boost/graph/adjacency_list.hpp>  // for source, target, etc
#include <boost/graph/properties.hpp>   // for default_color_type, etc
#include <boost/graph/topological_sort.hpp>  // for topological_sort
#include <boost/pending/property.hpp>   // for property
#include <boost/graph/graph_selectors.hpp>  // for directedS
#include <boost/graph/graph_traits.hpp>  // for graph_traits, etc

#include <iterator>                     // for back_insert_iterator, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc

namespace dealii
{
  namespace DoFRenumbering
  {
#ifndef DOXYGEN
    namespace boost
    {
      namespace boosttypes
      {
        using namespace ::boost;
        using namespace std;

        typedef adjacency_list<vecS, vecS, directedS,
        property<vertex_color_t, default_color_type,
        property<vertex_degree_t, int> > > Graph;
        typedef graph_traits<Graph>::vertex_descriptor Vertex;
        typedef graph_traits<Graph>::vertices_size_type size_type;

        typedef std::pair<size_type, size_type> Pair;
      }

      namespace internal
      {
        template <class DH>
        void
        downstream_graph (boosttypes::Graph &graph,
            std::vector<std::vector<bool> > &incoming_face,
            const DH &dof_handler,
            const Point<DH::space_dimension> &direction)
        {
          const unsigned int fe_degree = dof_handler.get_fe ().degree;
          QGauss<DH::space_dimension - 1> face_quadrature (fe_degree+1);
          const UpdateFlags face_update_flags = update_quadrature_points
          | update_normal_vectors;
          FEFaceValues<DH::space_dimension> fe_values_face (
              dof_handler.get_fe (), face_quadrature, face_update_flags);
          typename DH::active_cell_iterator cell = dof_handler.begin_active ();
          typename DH::active_cell_iterator endc = dof_handler.end ();
          for (; cell != endc; ++cell)
          {
            for (unsigned int face_no = 0;
                face_no < GeometryInfo<DH::space_dimension>::faces_per_cell;
                ++face_no)
            {
              // which assemblies share the previous edge
              typename DoFHandler<DH::space_dimension>::face_iterator face =
              cell->face (face_no);
              typename DoFHandler<DH::space_dimension>::cell_iterator neighbor =
              cell->neighbor (face_no);
              // Initializing the fe_values for the face and calling the normals
              fe_values_face.reinit (cell, face_no);
              const std::vector<Tensor<1,DH::space_dimension> > &normals =
              fe_values_face.get_all_normal_vectors ();
              // We only use one point because we assume straight faces
              unsigned int q_point = 0;
              // if the normal is in the same direction as the downstream,
              // we submit the entries to the boost graph
              if (direction * normals[q_point] > 0)
              {
                if (face->at_boundary ())
                continue;
                if (face->has_children ())
                {
                  for (unsigned int subface_no = 0;
                      subface_no < face->number_of_children (); ++subface_no)
                  {
                    add_edge (cell->user_index (),
                        cell->neighbor_child_on_subface (face_no,
                            subface_no)->user_index (), graph);
                  }
                }
                else
                {
                  add_edge (cell->user_index (),
                      cell->neighbor (face_no)->user_index (), graph);
                }
              }
              // if the normal is in oposite direction respet to the
              // downstream (or perpendicular) we say that it is an icoming face
              else
              {
                incoming_face[cell->user_index()][face_no] = true;
              }
            } //face_no
          } //cell
        } //downstream_graph
      }
    }
#endif

    template <class DH>
    void
    compute_downstream_cells (
        std::vector<typename DH::active_cell_iterator> & ordered_cells,
        std::vector<std::vector<bool> > & incoming_face,
        const DH &dof_handler,
        const Point<DH::space_dimension> &direction)
    {
      // we get the number of cells and resize the containers
#if DEAL_II_VERSION_GTE(8,4,0)
      unsigned int n_cells = dof_handler.get_triangulation ().n_active_cells ();
#else
      unsigned int n_cells = dof_handler.get_tria ().n_active_cells ();
#endif
      ordered_cells.resize (n_cells);
      incoming_face.resize (n_cells,
          std::vector<bool>(GeometryInfo<DH::space_dimension>::faces_per_cell, false));

      // initialize and fill the boost graph for the connections
      boost::boosttypes::Graph graph (n_cells);
      boost::internal::downstream_graph (graph, incoming_face, dof_handler, direction);

      // Then we reorder the graph with a topologycal sort
      typedef std::vector<boost::boosttypes::Vertex> container;
      container c;
      ::boost::topological_sort (graph, std::back_inserter (c));

      // here we should invert the ordering to pass to ordered_cells
      std::vector<unsigned int> invert_perm (n_cells);
      container::reverse_iterator it = c.rbegin ();
      for (unsigned int i = 0; i < n_cells; ++i, ++it)
      {
        invert_perm[*it] = i;
      }

      // store the ordered cells
      typename DH::active_cell_iterator cell = dof_handler.begin_active ();
      for (unsigned int i = 0; i < ordered_cells.size (); ++i, ++cell)
      {
        ordered_cells[invert_perm[i]] = cell;
      }
    }

    template <class DH>
    void
    compute_downstream_new (std::vector<types::global_dof_index> &new_indices,
        std::vector<types::global_dof_index> &reverse,
        const DH &dof_handler,
        const Point<DH::space_dimension> &direction)
    {
      // we get the number of cells and initialize the containers
#if DEAL_II_VERSION_GTE(8,4,0)
      unsigned int n_cells = dof_handler.get_triangulation ().n_active_cells ();
#else
      unsigned int n_cells = dof_handler.get_tria ().n_active_cells ();
#endif
      std::vector<typename DH::active_cell_iterator> ordered_cells (n_cells);
      std::vector<std::vector<bool> > dummy_incoming_face (n_cells,
          std::vector<bool>(GeometryInfo<DH::space_dimension>::faces_per_cell, false));

      // compute the ordered cells
      compute_downstream_cells (ordered_cells, dummy_incoming_face, dof_handler, direction);

      // this should be like this to move information from ordered cells
      // to new_indices and reverse, in order to be used with the sparse matrix
      compute_cell_wise (new_indices, reverse, dof_handler, ordered_cells);
    }

// Explicit code generation

#ifndef DOXYGEN
    template void
    compute_downstream_cells<DoFHandler<1> > (
        std::vector<typename DoFHandler<1>::active_cell_iterator> & ordered_cells,
        std::vector<std::vector<bool> > & incoming_face,
        const DoFHandler<1> & dof_handler,
        const Point<DoFHandler<1>::space_dimension> & direction);

    template void
    compute_downstream_cells<DoFHandler<2> > (
        std::vector<typename DoFHandler<2>::active_cell_iterator> & ordered_cells,
        std::vector<std::vector<bool> > & incoming_face,
        const DoFHandler<2> & dof_handler,
        const Point<DoFHandler<2>::space_dimension> & direction);

    template void
    compute_downstream_cells<DoFHandler<3> > (
        std::vector<typename DoFHandler<3>::active_cell_iterator> & ordered_cells,
        std::vector<std::vector<bool> > & incoming_face,
        const DoFHandler<3> & dof_handler,
        const Point<DoFHandler<3>::space_dimension> & direction);

    template void
    compute_downstream_new<DoFHandler<1> > (
        std::vector<types::global_dof_index> & new_indices,
        std::vector<types::global_dof_index> & reverse,
        const DoFHandler<1> & dof_handler,
        const Point<DoFHandler<1>::space_dimension> & direction);

    template void
    compute_downstream_new<DoFHandler<2> > (
        std::vector<types::global_dof_index> & new_indices,
        std::vector<types::global_dof_index> & reverse,
        const DoFHandler<2> & dof_handler,
        const Point<DoFHandler<2>::space_dimension> & direction);

    template void
    compute_downstream_new<DoFHandler<3> > (
        std::vector<types::global_dof_index> & new_indices,
        std::vector<types::global_dof_index> & reverse,
        const DoFHandler<3> & dof_handler,
        const Point<DoFHandler<3>::space_dimension> & direction);
#endif
  }
}
