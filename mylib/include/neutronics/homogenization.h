/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/homogenization.h
 * @brief  Homogenization class template declarations
 */

#ifndef FOREST_HOMOGENIZATION_H
#define FOREST_HOMOGENIZATION_H

#include "input/input_fwd.h"            // For class Input&
#include "input/input_mat.h"            // for InputMat::XsSingle
#include "neutronics/state_fwd.h"       // For class State&
#include "geometry/prob_geom_fwd.h"     // For class ProbGeom&

#include "deal.II/dofs/dof_handler.h"   // for DoFHandler
#include "deal.II/dofs/dof_accessor.h"  // for DoFHandler (necessary for 8.5)

#include <string>                       // for string
#include <vector>                       // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class ExtractBV; }
namespace dealii { template <int dim, int spacedim> class FiniteElement; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @brief Homogenization
   *
   * @ingroup ForestNeutronics
   *
   * @details The homogenization in space and condensation in energy
   * is performed by using the following definitions
   *
   * ## Homogenization and condensation of transport cross sections
   *
   * To understand the process of homogenization in space and condensation in
   * energy we can show the following two meshes
   *
   * \verbatim
   *            r1 r2 r3 r4
   *
   *            R1   R2  R3
   *           _____________
   * g1   G1   |  |  :  |  |
   *           _____________
   * g2        |  |  :  |  |
   *      G2   ----..-..----
   * g3        |  |  :  |  |
   *           _____________
   * g4   G3   |  |  :  |  |
   *           _____________
   * \endverbatim
   *
   * where {r1, r2, r3, r4} are the regions composing the fine mesh for the
   * space, and {g1, g2, g3, g4} is the fine mesh for the energy. The
   * coarsening is done by defining coarse regions {R1, R2, R3} and
   * {G1, G2, G3} in such a way that the coarse mesh match the fine mesh,
   * where the coarse regions contains at least one fine region (for example
   * R1 = {r1}), but possibly more (for example G2 = {g2, g3} ).
   *
   * Then we define the flux in a region \f$ r\f$ as the absolute value of
   * the flux for this region  (in space and energy)
   *
   * \f{align*}{
   * \Psi_{r,g} = \int_{V_{r}} \Psi_{g}(\mathbf{x}) \ \text{d}\mathbf{x}
   * \f}
   *
   * and the flux in the macro region (in space and energy)
   *
   * \f{align*}{
   * \Psi_{\mathbf{r},\mathbf{g}}
   * = \sum_{g \in \mathbf{g}} \sum_{r \in \mathbf{r}} \Psi_{r,g}
   * \f}
   *
   * Thus we use these quantities to perform a flux weighting of the cross
   * sections in the different regions (in space and energy) as follows
   *
   * \f{align*}{
   * \Sigma^{T}_{\mathbf{r},\mathbf{g}}
   * & = \frac{\sum_{g \in \mathbf{g}} \sum_{r \in \mathbf{r}}
   * \Sigma^{T}_{r,g} \Psi_{r,g} }
   * {\Psi_{\mathbf{r},\mathbf{g}} }
   * \f}
   *
   * \f{align*}{
   * \Sigma^{F}_{\mathbf{r},\mathbf{g}}
   * & = \frac{\sum_{g \in \mathbf{g}} \sum_{r \in \mathbf{r}}
   * \Sigma^{F}_{r,g} \Psi_{r,g} }{\Psi_{\mathbf{r},\mathbf{g}} }
   * \f}
   *
   * \f{align*}{
   * \nu \Sigma^{F}_{\mathbf{r},\mathbf{g}}
   * & = \frac{\sum_{g \in \mathbf{g}} \sum_{r \in \mathbf{r}}
   * \nu \Sigma^{F}_{r,g} \Psi_{r,g} }{\Psi_{\mathbf{r},\mathbf{g}} }
   * \f}
   *
   * \f{align*}{
   * \Sigma^{S}_{\mathbf{r},\mathbf{g} \to \mathbf{g}'}
   * & = \frac{\sum_{g \in \mathbf{g}} \sum_{g ' \in \mathbf{g} '}
   * \sum_{r \in \mathbf{r}} \Sigma^S_{r,g \to g '} \Psi_{r,g} }
   * {\Psi_{\mathbf{r},\mathbf{g}} }
   * \f}
   *
   * \f{align*}{
   * \chi_{\mathbf{r},\mathbf{g}}
   * & = \frac{\sum_{r \in \mathbf{r}} \sum_{g \in \mathbf{g}}
   * \chi_{r,g} \sum_{g ' = 1}^G \nu \Sigma^F_{r,g '} \Psi_{r,g '} }
   * {\sum_{r \in\mathbf{r}}  \sum_{g '' = 1}^G
   * \chi_{r,g ''} \sum_{g ' = 1}^G \nu \Sigma^{F}_{r,g '} \Psi_{r,g '} }
   * = \frac{ \sum_{r \in \mathbf{r}} f_{r}
   * \sum_{g \in \mathbf{g}} \chi_{r,g} }
   * { \sum_{r \in\mathbf{r}} f_{r}
   * \underbrace{ \sum_{g '' = 1}^G \chi_{r,g ''} }_{=1} } \\
   * & = \sum_{r \in \mathbf{r}} \frac{ f_{r} }
   * { \sum_{r \in\mathbf{r}} f_{r}  } \sum_{g \in \mathbf{g}}  \chi_{r,g}
   * = \sum_{r \in \mathbf{r}} \frac{f_{r}}{\sum_{r \in\mathbf{r}} f_{r} } \chi_{r,\mathbf{g}}
   * \f}
   *
   * \f{align*}{
   * \Sigma^{tr}_{\mathbf{r},g}
   * & = \frac{\sum_{r \in \mathbf{r}} \Sigma^{tr}_{r,g} \Psi_{r,g} }
   * {\sum_{r \in \mathbf{r}} \Psi_{r,g} }
   * \f}
   *
   * ## Homogenization and condensation of diffusion coefficient
   *
   * The diffusion coefficient is usually homogenized, which is defined as
   * one third of the inverse of the transport cross section, should be
   * considered a bit in detail.
   *
   * It is calculated by first homogenizing in space the transport cross section,
   *
   * \f{align*}{
   * D_{\mathbf{r},g} = \frac{1}{3 \Sigma^{tr}_{\mathbf{r},g}}
   * \f}
   *
   * and then the standard formula for the energy condensation
   *
   * \f{align*}{
   * D_{\mathbf{r},\mathbf{g}}
   * = \frac{\sum_{g \in \mathbf{g}} D_{\mathbf{r},g}
   * \Psi_{\mathbf{r},g}}{\Psi_{\mathbf{r},\mathbf{g}}}
   * = \frac{\sum_{g \in \mathbf{g}} \frac{1}{3 \Sigma^{tr}_{\mathbf{r},g}}
   * \Psi_{\mathbf{r},g}}{\Psi_{\mathbf{r},\mathbf{g}}}
   * \f}
   *
   * ## Absorbtion cross section (to be used with the diffusion equation)
   *
   * The absorbtion cross section \f$ \Sigma^A_{\mathbf{g}} \f$ is defined as
   *
   * \f{align*}{
   * \Sigma^A_{{g}}
   * \equiv \Sigma^T_{g}
   * - \sum_{\substack{g' = 1}}^{G} \Sigma^S_{g \to g'}
   * \f}
   *
   * ## Removal cross section (to be used with the diffusion equation)
   *
   * The removal cross section \f$ \Sigma^R_{\mathbf{g}} \f$ is defined as
   *
   * \f{align*}{
   * \Sigma^R_{{g}}
   * \equiv \Sigma^A_{g}
   * + \sum_{\substack{g' = 1 \\ g' \neq g}}^{G} \Sigma^S_{g \to g'}
   * \f}
   *
   * or alternatively, using the definition of the absortion cross section, as
   *
   * \f{align*}{
   * \Sigma^R_{{g}}
   * & \equiv \Sigma^A_{g}
   * + \sum_{\substack{g' = 1 \\ g' \neq g}}^{G} \Sigma^S_{g \to g'} \\
   * & = \Sigma^T_{g}
   * - \sum_{\substack{g' = 1}}^{G} \Sigma^S_{g \to g'}
   * + \sum_{\substack{g' = 1 \\ g' \neq g}}^{G} \Sigma^S_{g \to g'} \\
   * & = \Sigma^T_{g} - \Sigma^S_{g \to g}
   * \f}
   *
   * If axial buckling is specified the removal cross section writes as
   *
   * \f{align*}{
   * \Sigma^R_{g}
   * \equiv D_{g} B_z^2 + \Sigma^A_{g}
   * + \sum_{\substack{g' = 1 \\ g' \neq g}}^{G} \Sigma^S_{g \to g'}
   * \f}
   *
   *
   */
  template <int dim>
  class Homogenization
  {
  public:

    /** @brief Constructor */
    Homogenization (const State<dim> &state/*,
     const ExtractBV<dim> & extract_bv*/);

    /** @brief Calculate flux weight per cell */
    void
    cell_weights ();

    /** @brief Homogenize at pin level */
    void
    init_xs (XS_single & dst,
             const unsigned int n_groups);

    /** @brief Homogenize at pin level */
    void
    add_cell_contribution (XS_single & dst,
                           const XS_single & src,
                           const std::vector<double> & weight);

    /** @brief Homogenize at pin level */
    void
    normalize (XS_single & dst,
               const std::vector<double> &  weight);

    /** @brief Homogenize at pin level */
    void
    pin_homogenization ();

    /** @brief Homogenize at assembly level */
    void
    assembly_homogenization ();

  private:

    /** @brief Short name for the face iterator. */
    typedef typename DoFHandler<dim, dim>::face_iterator f_iter;

    /** @brief Short name for the cell iterator. */
    typedef typename DoFHandler<dim, dim>::active_cell_iterator c_iter;

    const State<dim> & mp_state;
    /*const ExtractBV<dim> & mp_extract_bv;*/
    const Input & data;
    const ProbGeom<dim> & geom;
    const FiniteElement<dim, dim> & fe;
    const DoFHandler<dim, dim> & dof_handler;

    std::vector<std::vector<double> > cell_w;
    std::vector<std::vector<double> > pin_w;
    std::vector<std::vector<double> > ass_w;

  public:
    std::vector<XS_single> pin_xs;
    std::vector<XS_single> ass_xs;
  private:

    /** @brief shortcut for std::vector<double> */
    typedef std::vector<double> stdvd1;
    /** @brief recursive shortcut for std::vector<std_vd1> */
    typedef std::vector<stdvd1> stdvd2;
    /** @brief recursive shortcut for std::vector<std_vd2> */
    typedef std::vector<stdvd2> stdvd3;
    /** @brief recursive shortcut for std::vector<std_vd3> */
    typedef std::vector<stdvd3> stdvd4;

    /** @brief Print the homogenized cross sections */
    /*void
     print_homogenized_xs (const std::string &filename,
     const std::vector<InputMat::XS_single> & xs,
     const stdvd3 & phi0,
     const stdvd3 & phi2,
     const stdvd3 & current);*/
  };

} /* end of namespace ForestNeutronics */

#endif /* FOREST_HOMOGENIZATION_H */
