// ---------------------------------------------------------------------
// @f$Id@f$
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef path_to_intergridmap_sebas_h
#define path_to_intergridmap_sebas_h

#include "deal.II/base/exceptions.h"    // for DeclException0, etc
#include "deal.II/base/smartpointer.h"  // for SmartPointer
#include "deal.II/base/subscriptor.h"   // for Subscriptor

#include <ostream>                      // for operator<<, basic_ostream
#include <vector>
#include <cstddef>

#ifndef DOXYGEN
namespace dealii { class ConstraintMatrix; }
#endif

namespace dealii
{
  /**
   * @brief InterGridMapSebas
   */
  template <class GridClass>
  class InterGridMapSebas : public Subscriptor
  {
  public:

    /** @todo document me */
    typedef typename GridClass::cell_iterator cell_iterator;

    /** @todo document me */
    InterGridMapSebas ();

    /** @todo document me */
    void
    make_mapping (const GridClass &source_grid,
                  const GridClass &destination_grid);

    /** @todo document me */
    cell_iterator
    operator [] (const cell_iterator &source_cell) const;

    /** @todo document me */
    void
    clear ();

    /** @todo document me */
    const GridClass &
    get_source_grid () const;

    /** @todo document me */
    const GridClass &
    get_destination_grid () const;

    /** @todo document me */
    std::size_t
    memory_consumption () const;

    /** @todo document me */
    DeclException1 (ExcInvalidKey,
        cell_iterator,
        << "The iterator " << arg1 << " is not valid as key for "
        << "this map.")
    ;

    /** @todo document me */
    DeclException0 (ExcIncompatibleGrids)
    ;

  private:
    std::vector<std::vector<cell_iterator> > mapping;

    SmartPointer<const GridClass, InterGridMapSebas<GridClass> > source_grid;

    SmartPointer<const GridClass, InterGridMapSebas<GridClass> > destination_grid;

    void
    set_mapping (const cell_iterator &src_cell, const cell_iterator &dst_cell);

    void
    set_entries_to_cell (const cell_iterator &src_cell,
                         const cell_iterator &dst_cell);
  };

  namespace VectorTools
  {
    /**
     * @brief interpolate_to_different_mesh
     */
    template <int dim, int spacedim, template <int, int> class DH, class VECTOR>

    /** @todo document me */
    void
    interpolate_to_different_mesh (
        const InterGridMapSebas<DH<dim, spacedim> > &intergridmap,
        const VECTOR &u1, const ConstraintMatrix &constraints, VECTOR &u2);
  }
}
#endif
