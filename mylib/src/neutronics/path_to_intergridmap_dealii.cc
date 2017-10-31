/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/path_to_intergridmap_dealii.cc
 * @brief  Implementation of grid generator functions
 */

//---------------------------------------------------------------------------
//    $Id: intergrid_map.cc 23709 2011-05-17 04:34:08Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,
//                  2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "neutronics/path_to_intergridmap_dealii.h"

#include "deal.II/base/config.h"
#include "deal.II/base/memory_consumption.h"
#include "deal.II/base/smartpointer.h"
#include "deal.II/grid/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/dofs/dof_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/intergrid_map.h"
#include "deal.II/dofs/dof_accessor.h"

namespace dealii
{

  namespace
  {
// helper function to aquire the number of levels within a grid
    template <class GridClass>
    unsigned int
    get_n_levels (const GridClass &grid)
    {
      // all objects we deal with are able
      // to deliver a pointer to the
      // underlying triangulation.
      //
      // for the triangulation as GridClass
      // of this object, there is a
      // specialization of this function
#if DEAL_II_VERSION_GTE(8,4,0)
      return grid.get_triangulation ().n_levels ();
#else
      return grid.get_tria ().n_levels ();
#endif
    }

// specialization for grid==tria
    template <int dim>
    unsigned int
    get_n_levels (const Triangulation<dim, dim> &grid)
    {
      // if GridClass==Triangulation, then
      // we can ask directly.
      return grid.n_levels ();
    }
  }

  template <class GridClass>
  InterGridMapSebas<GridClass>::InterGridMapSebas ()
      : source_grid (0, typeid(*this).name ()),
        destination_grid (0, typeid(*this).name ())
  {
  }

  template <class GridClass>
  void
  InterGridMapSebas<GridClass>::make_mapping (const GridClass &source_grid,
                                              const GridClass &destination_grid)
  {
    // first delete all contents
    clear ();

    // next store pointers to grids
    this->source_grid = &source_grid;
    this->destination_grid = &destination_grid;

    // then set up the containers from
    // scratch and fill them with end-iterators
    const unsigned int n_levels = get_n_levels (source_grid);
    mapping.resize (n_levels);
    for (unsigned int level = 0; level < n_levels; ++level)
    {
      // first find out about the highest
      // index used on this level. We could
      // in principle ask the triangulation
      // about this, but we would have to
      // know the underlying data structure
      // for this and we would like to
      // avoid such knowledge here
      unsigned int n_cells = 0;
      cell_iterator cell = source_grid.begin (level), endc = source_grid.end (
          level);
      for (; cell != endc; ++cell)
      {
        if (static_cast<unsigned int> (cell->index ()) > n_cells)
        {
          n_cells = cell->index ();
        }
      }

      // note: n_cells is now the largest
      // zero-based index, but we need the
      // number of cells, which is one larger
      mapping[level].resize (n_cells + 1, destination_grid.end ());
    };

    // now make up the mapping
    // loop over all cells and set the user
    // pointers as well as the contents of
    // the two arrays. note that the function
    // takes a *reference* to the int and
    // this may change it
    cell_iterator src_cell = source_grid.begin (0), dst_cell =
        destination_grid.begin (0), endc = source_grid.end (0);
    for (; src_cell != endc; ++src_cell, ++dst_cell)
    {
      set_mapping (src_cell, dst_cell);
    }

    // little assertion that the two grids
    // are indeed related:
    Assert(dst_cell == destination_grid.end(0), ExcIncompatibleGrids ());
  }

  template <class GridClass>
  void
  InterGridMapSebas<GridClass>::set_mapping (const cell_iterator &src_cell,
                                             const cell_iterator &dst_cell)
  {
    // first set the map for this cell
    mapping[src_cell->level ()][src_cell->index ()] = dst_cell;

    // if both cells have children, we may
    // recurse further into the hierarchy
    if (src_cell->has_children () && dst_cell->has_children ())
    {
      /* sebas modifications commenting the previous work
       * Assert(src_cell->n_children()==
       * GeometryInfo<GridClass::dimension>::max_children_per_cell,
       * ExcNotImplemented());
       * Assert(dst_cell->n_children()==
       * GeometryInfo<GridClass::dimension>::max_children_per_cell,
       * ExcNotImplemented());
       */
      // now sebas modifications
      Assert(
          src_cell->n_children() <= GeometryInfo<GridClass::dimension>::max_children_per_cell,
          ExcNotImplemented());
      Assert(
          dst_cell->n_children() <= GeometryInfo<GridClass::dimension>::max_children_per_cell,
          ExcNotImplemented());
      // until here

      Assert(src_cell->refinement_case()==dst_cell->refinement_case(),
          ExcNotImplemented());

      // now my modifications
      // for (unsigned int c=0; c<GeometryInfo<GridClass::dimension>::max_children_per_cell; ++c)
      //   set_mapping (src_cell->child(c),
      //                dst_cell->child(c));
      for (unsigned int c = 0; c < src_cell->n_children (); ++c)
      {
        set_mapping (src_cell->child (c), dst_cell->child (c));
      }
    }
    else if (src_cell->has_children () && !dst_cell->has_children ())
    {
      // src grid is more refined here.
      // set entries for all children
      // of this cell to the one
      // dst_cell
      for (unsigned int c = 0; c < src_cell->n_children (); ++c)
      {
        set_entries_to_cell (src_cell->child (c), dst_cell);
      }
    }
    // else (no cell is refined or
    // dst_cell is refined): no pointers
    // to be set
  }

  template <class GridClass>
  void
  InterGridMapSebas<GridClass>::set_entries_to_cell (const cell_iterator &src_cell,
                                                     const cell_iterator &dst_cell)
  {
    // first set the map for this cell
    mapping[src_cell->level ()][src_cell->index ()] = dst_cell;

    // then do so for the children as well
    // if there are any
    if (src_cell->has_children ())
    {
      for (unsigned int c = 0; c < src_cell->n_children (); ++c)
      {
        set_entries_to_cell (src_cell->child (c), dst_cell);
      }
    }
  }

  template <class GridClass>
  typename InterGridMapSebas<GridClass>::cell_iterator
  InterGridMapSebas<GridClass>::operator [] (
      const cell_iterator &source_cell) const
  {
    Assert(source_cell.state() == IteratorState::valid,
        ExcInvalidKey (source_cell));
    Assert(source_cell->level() <= static_cast<int>(mapping.size()),
        ExcInvalidKey (source_cell));
    Assert(
        source_cell->index() <= static_cast<int>(mapping[source_cell->level()].size()),
        ExcInvalidKey (source_cell));

    return mapping[source_cell->level ()][source_cell->index ()];
  }

  template <class GridClass>
  void
  InterGridMapSebas<GridClass>::clear ()
  {
    mapping.clear ();
    source_grid = 0;
    destination_grid = 0;
  }

  template <class GridClass>
  const GridClass &
  InterGridMapSebas<GridClass>::get_source_grid () const
  {
    return *source_grid;
  }

  template <class GridClass>
  const GridClass &
  InterGridMapSebas<GridClass>::get_destination_grid () const
  {
    return *destination_grid;
  }

  template <class GridClass>
  std::size_t
  InterGridMapSebas<GridClass>::memory_consumption () const
  {
    return (MemoryConsumption::memory_consumption (mapping)
        + MemoryConsumption::memory_consumption (source_grid)
        + MemoryConsumption::memory_consumption (destination_grid));
  }

// explicit instantiations
  template class InterGridMapSebas<DoFHandler<1, 1> > ;
  template class InterGridMapSebas<DoFHandler<2, 2> > ;
  template class InterGridMapSebas<DoFHandler<3, 3> > ;

  namespace VectorTools
  {

    template <int dim, int spacedim, template <int, int> class DH, class Vector>
    void
    interpolate_to_different_mesh (const InterGridMapSebas<DH<dim, spacedim> > &intergridmap,
                                   const Vector &u1,
                                   const ConstraintMatrix &constraints,
                                   Vector &u2)
    {
      const DH<dim, spacedim> &dof1 = intergridmap.get_source_grid ();
#ifdef DEBUG
      const DH<dim, spacedim> &dof2 = intergridmap.get_destination_grid ();
#endif

      Assert(u1.size()==dof1.n_dofs(),
          ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
      Assert(u2.size()==dof2.n_dofs(),
          ExcDimensionMismatch(u2.size(),dof2.n_dofs()));

      Vector cache;

      // Looping over the finest common
      // mesh, this means that source and
      // destination cells have to be on the
      // same level and at least one has to
      // be active.
      //
      // Therefor, loop over all cells
      // (active and inactive) of the source
      // grid ..
      typename DH<dim, spacedim>::cell_iterator cell1 = dof1.begin ();
      const typename DH<dim, spacedim>::cell_iterator endc1 = dof1.end ();

      for (; cell1 != endc1; ++cell1)
      {
        const typename DH<dim, spacedim>::cell_iterator cell2 =
            intergridmap[cell1];

        // .. and skip if source and destination
        // cells are not on the same level ..
        if (cell1->level () != cell2->level ())
        {
          continue;
        }
        // .. or none of them is active.
        if (!cell1->active () && !cell2->active ())
        {
          continue;
        }

        Assert(cell1->get_fe().get_name() == cell2->get_fe().get_name(),
            ExcMessage ("Source and destination cells need to use the same finite element"));

        cache.reinit (cell1->get_fe ().dofs_per_cell);

        // Get and set the corresponding
        // dof_values by interpolation.
        cell1->get_interpolated_dof_values (u1, cache);
        cell2->set_dof_values_by_interpolation (cache, u2);
      }

      // Apply hanging node constraints.
      constraints.distribute (u2);
    }

#ifndef DOXYGEN
// explicit instantiations

    template
    void
    interpolate_to_different_mesh (const InterGridMapSebas<DoFHandler<1, 1> > &,
                                   const Vector<double> &,
                                   const ConstraintMatrix &,
                                   Vector<double> &);

    template
    void
    interpolate_to_different_mesh (const InterGridMapSebas<DoFHandler<2, 2> > &,
                                   const Vector<double> &,
                                   const ConstraintMatrix &,
                                   Vector<double> &);

    template
    void
    interpolate_to_different_mesh (const InterGridMapSebas<DoFHandler<3, 3> > &,
                                   const Vector<double> &,
                                   const ConstraintMatrix &,
                                   Vector<double> &);
#endif
  }

}
