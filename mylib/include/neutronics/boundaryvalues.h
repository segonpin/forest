/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/boundaryvalues.h
 * @brief  BoundaryValues class template declarations
 */

#ifndef FOREST_BOUNDARYVALUES_H
#define FOREST_BOUNDARYVALUES_H

/**
 * @warning We need "state.h" included. If we remove it this will not compile
 * for me (sebas) due to undefined behavior of one of the template parameters
 * for the multimap (c_iter)  */
#include "neutronics/state.h"
#include "angle/quadrature_base_fwd.h"
#include "neutronics/state_fwd.h"

#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/dofs/dof_accessor.h"  // for DoFHandler (necessary for 8.5)

#include <vector>                       // for vector
#include <map>                          // for multimap, map
#include <memory>                       // for shared_ptr

#ifndef DOXYGEN
namespace dealii { template <typename number> class Vector; }
namespace dealii { template <typename number> class BlockVector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**  
   @ingroup ForestNeutronics
   Boundary conditions for the system.

   After the discretization of the energy and the angle (in this case with
   discrete ordinates), the streaming operator over a domain \f$ T \f$ has the
   following form

   \f{align*}{
   & \Omega_m \cdot \nabla \Psi_{m,g}(\mathbf{r})
   \quad \forall m, g \quad \mathbf{r} \in T .
   \f}

   And when the finite element method is used, using integration by parts we get

   \f{align*}{
   & \int_T \left( \Omega_m \cdot \nabla \right)
   \Psi_{m,g}(\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r}) \\
   & = \int_T  \Psi_{m,g}(\mathbf{r})
   \left( \Omega_m \cdot \nabla \right) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   - \int_{\partial T}\left( \Omega_m \cdot \vec{n}(\mathbf{r}) \right)
   \Psi_{m,g}(\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r})  \\
   & = \int_T \Psi_{m,g}(\mathbf{r})
   \left( \Omega_m \cdot \nabla \right) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   - \int_{\partial T_{\text{out}}}
   \left( \Omega_m \cdot \vec{n}(\mathbf{r}) \right)
   \Psi_{m,g}(\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   - \int_{\partial T_{\text{in}}}
   \left( \Omega_m \cdot \vec{n}(\mathbf{r}) \right)
   \Psi_{m,g}(\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   \f}

   Then, the boundary conditions are imposed by substituting the value of
   \f$ \Psi_{m,g}(\mathbf{r}) \f$ at the incoming boundaries, by its value,
   and moving this term to the right hand side (when the value is known), or
   generating coupling between directions when this value is related to an
   outgoing flux for a different direction (reflective boundary conditions).

   ## Incoming (or Dirichlet) Boundary Conditions.

   Incoming  boundary conditions can be written as

   \f{align*}{
   & \Psi_{m,g}(\mathbf{r}) = \Psi_{m,g}^{\text{in}}(\mathbf{r}),
   \quad n(\mathbf{r})\cdot \Omega_m < 0 \quad \forall m, g
   \quad \forall \mathbf{r} \in \Gamma_{D}
   \f}

   being \f$ \Gamma_{D} \f$ the incoming part of the boundary (Dirichlet).
   In matrix form, this can be used to solve the system

   \f{align*}{
   L_{0,g} \Psi_{g} = Q_{g} - L_{D,g} \Psi_{g}^{\text{in}}
   \f}

   where \f$ L_{0,g} \f$ is a block diagonal matrix, where the blocks in the
   diagonal are constructed as the discretized version of the **bilinear** form

   \f{align*}{
   \mathcal{L}_{0,m,g} (\Psi_{m,g},v) = \int_T \Psi_{m,g}(\mathbf{r})
   \left( \Omega_m \cdot \nabla \right) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   + \int_{\partial T_{\text{out}}}
   \left( \Omega_m \cdot \vec{n}(\mathbf{r}) \right)
   \Psi_{m,g}(\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   \quad \forall m
   \f}

   and the right hand side is the discretized version of the **linear** forms

   \f{align*}{
   Q_{m,g}(v) - \mathcal{L}_{D,m,g} (v)
   = \int_T Q_{m,g}(\mathbf{r})
   v(\mathbf{r}) \ \text{d}(\mathbf{r})
   + \int_{\partial T_{\text{in}}}
   \left( \Omega_m \cdot \vec{n}(\mathbf{r}) \right)
   \Psi_{m,g}^{\text{in}} (\mathbf{r}) v(\mathbf{r}) \ \text{d}(\mathbf{r})
   \quad \forall m
   \f}

   ## Reflective Boundary Conditions.

   Reflective boundary conditions can be written as

   \f{align*}{
   & \Psi_{m,g}(\mathbf{r}) = \Psi_{m,g}^{\text{in}}(\mathbf{r})
   := \beta_{m',g} \Psi_{m',g}^{\text{out}}(\mathbf{r}) ,
   \quad n(\mathbf{r})\cdot \Omega_m < 0 \quad \forall m,g
   \quad \forall \mathbf{r} \in \Gamma_{R}
   \f}

   being \f$ \Gamma_{R}\f$ the boundary with reflecting boundary conditions.
   It is important to notice that the index \f$ m' \f$ representing the
   outgoing direction that will enter as an incoming direction for \f$ m \f$
   is calculated by

   \f{align*}{
   \Omega_{m'} := \Omega_m - 2 \vec{n}(\mathcal{r})
   (\Omega_m \cdot \vec{n}(\mathcal{r}))
   \f}

   and directly depends on the normal vector to the particular face, so it will
   be different for every face orientation.

   Here we drop the group index for simplicity.
   In order to clarify the coupling between the different directions, we can
   write the complete streaming operator containing all the directions as

   \f{align*}{
   L :=
   \left( \begin{array}{cccc}
   L_{0,m} & \ldots & \ldots & \\
   \ldots & L_{0,m} & \beta_{m} L_{R,m'} & \ldots \\
   \ldots & \beta_{n'} L_{R,n'} & L_{0,n} & \ldots \\
   & \ldots & \ldots  &  L_{0,M}
   \end{array} \right)
   = L_{0} + L_{R}
   \f}

   meaning that using one of the rows of the matrix L and multiplying for the
   vector of all the directions we get

   \f{align*}{
   L_{0,m} \Psi_{m} + \beta_{m'} L_{R,m'} \Psi_{m'}
   + \beta_{m''} L_{R,m''} \Psi_{m''} + \ldots
   \f}

   Here we notice that one direction can have more than one outgoing direction
   associated, with a maximum of one per each different oriented face of the
   boundary.

   We also notice that the boundary values, when lagged from a previous
   iteration, will be stored in a compressed form, meaning that only the values
   related to significant degrees of freedom are stored, representing them as
   \f$ \Psi_{S} \f$, and they can be obtained by applying the projection matrix
   \f$ P \f$ to the vector containing all the degrees of freedom, as follows

   \f{align*}{
   \Psi_S = P^{T}*\Psi , \quad  \iff \quad \Psi = P*\Psi_{S}
   \f}

   so the matrix \f$ L_R \f$ should be applied as follows

   \f{align*}{
   - L_{R}*\Psi := - L_{R}*P*\Psi_{S} .
   \f}

   @note The notation is inspired from @cite warsa2004krylov .

   @note It is important to notice here that the relations between the
   incoming and outgoing directions does not depends on the group,
   discretization. Nevertheless, the notation is maintained because the incoming
   flux, if provided, will most probably depend on each group, and the albedo
   factor applied to the reflective boundary conditions too.

   @todo do I need this?
   @code{.cpp}
   void get_albedos ();
   @endcode

   */
  template <int dim>
  class BoundaryValues
  {
  public:
    /**
     @brief Constructor using the information in @p state
     @param state
     @param ord
     */
    BoundaryValues (State<dim> & state,
                    std::shared_ptr<QuadratureBase<dim> > & ord);

  private:

    /** @brief Short name for the face iterator. */
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;

    /** @brief Short name for the cell iterator. */
    typedef typename DoFHandler<dim, dim>::active_cell_iterator c_iter;

    /** @brief Pointer to state object. */
    State<dim> & mp_state;

    /** @brief Pointer to quadrature object. */
    std::shared_ptr<QuadratureBase<dim> > & mp_ord;

    /**
     @brief For running over the dofs of the face and skipping the interiors.
     @details For the indices f,i representing the face and the local dof
     for the face, it returns the dofs for the cell.
     */
    std::vector<std::vector<unsigned int> > m_fe_sup_on_face;

    /**
     @brief Container for the incoming faces.
     This is intended to be used with an iterator for each direction @c i
     @code{.cpp}
     typename std::multimap<c_iter, unsigned int>::iterator it;
     for (it = m_i_faces[i].begin (); it != m_i_faces[i].end (); ++it)
     {
     c_iter cell = it->first;
     unsigned int face_no = it->second;
     }
     @endcode
     */
    std::vector<std::multimap<c_iter, unsigned int> > m_i_faces;

    /** @brief Container for the outgoing faces. See @ref m_i_faces . */
    std::vector<std::multimap<c_iter, unsigned int> > m_o_faces;

    /**
     @brief The matrix \f$ P \f$ extracting significant dofs for bcs.
     @details For a given direction @c dir , and a given dofs @c i,
     it gives me the significant index @c i_s to collect the dofs
     @code{.cpp}
     i_s = m_global_to_bc_data[dir][i]
     psi_s[dir][i_s] == psi[dir][i]
     @endcode
     @note The notation is a short version of m_global_to_significant_dofs
     */
    std::vector<std::map<unsigned int, unsigned int> > m_g2s_dofs;

    /** @brief Significant dofs for the incoming flux for every direction. */
    std::vector<unsigned int> m_dofs_in_bc;

  public:

    /**
     @brief apply the boundary conditions in @p psi_s
     for direction @p i_ord to the flux @p src_m
     @details For group @p g, @p psi_s contains the significant dofs for the
     boundary conditions, and these are applied to the vector
     @p src_m (\f$ Q_m \f$)

     \f{align*}{
     Q_{m,g} = Q_{m,g} - L_{R,m',g}*\Psi_{m',g} - L_{R,m'',g}*\Psi_{m'',g}
     - \ldots ,
     \f}

     but it does using the vector of the significant boundary values
     @p phi_s (\f$ \Psi_{S,g} \f$) as follows

     \f{align*}{
     Q_{m,g} = Q_{m,g} - L_{R,g}*P*\Psi_{S,g} .
     \f}

     @note The notation for \f$ -L_{R}*P*\psi_R \f$ is coming from the article
     @cite warsa2004krylov
     : "Krylov iterative methods applied to multidimensional Sn calculations in
     the presence of material discontinuities".
     @param src_m
     @param psi_s         The significant values of the boundary conditions.
     @param g             The energy group.
     @param i_dir         The incoming direction.
     @param f_bc_ind = 0  First index for the boundary condition.
     */
    void
    apply_reflective (Vector<double> & src_m,
                      const BlockVector<double> & psi_s,
                      unsigned int g,
                      unsigned int i_dir,
                      unsigned int f_bc_ind = 0);

    /**
     @brief Collect the outgoing flux for direction @p o_dir to incoming flux.
     @details Get the outgoing flux in direction @p o_dir and store it in
     @p psi_s for the corresponding incoming direction when it applies.
     The notation for this function would be the operator  \f$ P^{T} \f$, which
     applied to the vector with the degrees of freedoms for all the fluxes
     restrict the vector to only the significant degrees of freedom, obtaining

     \f{align*}{
     \Psi_{S,g} = P^{T}\Psi_{g} ,
     \f}

     Because it is applied for a vector with a particular outgoing direction,
     it actually the corresponding incoming direction for every type of face,
     and store this part of the significant degrees of freedom

     \f{align*}{
     \left \{ \begin{array}{l}
     \Psi_{S,m'} = P_{m'}^{T}\Psi^{\text{in}}_{m',g} :=
     P_{m'}^{T}\Psi^{\text{out}}_{m,g} , \\
     \Psi_{S,m''} = P_{m''}^{T}\Psi^{\text{in}}_{m'',g} :=
     P_{m''}^{T}\Psi^{\text{out}}_{m,g} , \\
     \quad \vdots
     \end{array} \right .
     \f}

     Where in this case the outgoing direction \f$ m \f$ has the incoming
     directions \f$ m', m'', \ldots \f$ .

     @note The notation for \f$ P^{T} \f$ is coming from the article
     @cite warsa2004krylov
     : "Krylov iterative methods applied to multidimensional
     Sn calculations in the presence of material discontinuities".
     @param psi_s Vector with the boundary conditions.
     @param psi_m
     @param o_dir  Outgoing direction of the flux that will be redistributed
     in the corresponding ingoing directions.
     @param f_bc_ind = 0. First index for the boundary condition.

     @todo  I am not sure about this line ...
     @code{.cpp}
     if (face->boundary_id () != f)
       continue;
     @endcode
     */
    void
    get_reflective (BlockVector<double> & psi_s,
                    const Vector<double> & psi_m,
                    unsigned int o_dir,
                    unsigned int f_bc_ind = 0);

    /**
     @brief Calculate memory consumption for this class.
     @return Return the memory consumption.
     */
    double
    memory_consumption () const
    {
      double memory_consumption = 0;
      memory_consumption += sizeof(m_i_faces) + sizeof(m_o_faces)
                            + sizeof(m_g2s_dofs)
                            + sizeof(m_dofs_in_bc);
      return memory_consumption;
    }

    /**
     @brief Significant dofs for the incoming flux for every direction.
     @return A @p n_angles vector with the significant dofs for each direction.
     */
    const std::vector<unsigned int> &
    get_block_sizes () const
    {
      return m_dofs_in_bc;
    }

    /**
     @brief Return the number of degrees of freedom for direction @p i_ord.
     @param i_ord
     @return significant dofs for direction @p i_ord.
     */
    unsigned int
    get_block_sizes (unsigned int i_ord) const
    {
      return m_dofs_in_bc[i_ord];
    }

  private:

    /**
     @brief Initializing the data. Necessary before using the class.
     @details Here we use the information in quadrature and geom to
     set which faces are incoming faces for a particular direction,
     and the dofs are saved to apply the bcs faster in the next iterations.
     */
    void
    init ();

    /**
     @brief Collect the albedo faces for direction i_ord.
     @details When a face is albedo, we store the cell iterator of this face,
     together with the index representing this face for the cell, in a vector,
     for a particular direction.
     @param i_ord
     */
    void
    set_incoming_faces (unsigned int i_ord);

    /**
     @brief Collect the albedo faces for direction i_ord.
     @details Similar to @ref set_incoming_faces. For the moment, it is used
     only for checking that the ingoing faces are the same number of the
     outgoing faces.
     @param i_ord
     */
    void
    set_outgoing_faces (unsigned int i_ord);

    /**
     @param i_ord
     */
    void
    set_global_to_bc_data (unsigned int i_ord);

  };

} // end of namespace Forest
#endif /* FOREST_BOUNDARYVALUES_H */
