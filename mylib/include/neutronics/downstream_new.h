/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/downstream_new.h
 * @brief  DofRenumbering utilities
 */

#ifndef FOREST_DOWNSTREAM_NEW_H
#define FOREST_DOWNSTREAM_NEW_H

#include "deal.II/base/types.h"
#include "deal.II/base/point.h"

#include <vector>

namespace dealii
{
  namespace DoFRenumbering
  {
    /** @brief compute_downstream_cells */
    template <class DH>
    void
    compute_downstream_cells (std::vector<typename DH::active_cell_iterator> & ordered_cells,
                              std::vector<std::vector<bool> > & incoming_face,
                              const DH &dof_handler,
                              const Point<DH::space_dimension> &direction);

    /** @brief compute_downstream_new */
    template <class DH>
    void
    compute_downstream_new (std::vector<types::global_dof_index> &new_indices,
                            std::vector<types::global_dof_index> &reverse,
                            const DH &dof,
                            const Point<DH::space_dimension> &direction);
  }
}

#endif /* FOREST_DOWNSTREAM_NEW_H */
