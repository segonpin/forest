/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/equation.h
 * @brief   Equation class template declarations
 */

#ifndef FOREST_EQUATION_H
#define FOREST_EQUATION_H

#include "algebra/sn_vector_fwd.h"

#include <vector>
#include <memory>

#ifndef DOXYGEN
namespace dealii { template <typename > class BlockVector; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class Equation
   * @brief Equation, base class
   * @ingroup ForestAlgebra
   *
   * @details Base class to construct the derived equations type. This class
   * specifies the contract to be satisfied by the derived class.
   *
   * The equation needs to know three different sizes, in order to generate
   * temporary vectors of the right size for the solvers. The three different
   * sizes are:
   *
   *   - Angular flux size (for recovering the full size).
   *   - Legendre moments + significant bdry dofs (for sg-gmres).
   *   - Fission density size (for eigen problem and for fission density).
   *
   */
  template <int dim>
  class Equation
  {
  public:

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_angular_sizes () const = 0;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_flux_bv_sizes () const = 0;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_fission_density_sizes () const = 0;

    /** @brief Pure virtual function defining the matrix vector product. */
    virtual void
    vmult (SnVector &dst,
           const SnVector &src) const = 0;

    /** @brief Pure virtual function defining the matrix vector product. */
    virtual void
    vmult (BlockVector<double> &dst,
           const BlockVector<double> &src) const = 0;

    /** @brief Pure virtual function. */
    virtual void
    Tvmult (SnVector & dst,
            const SnVector & src) const = 0;

    /** @brief Pure virtual function. */
    virtual void
    Tvmult (BlockVector<double> & dst,
            const BlockVector<double> & src) const = 0;

    /** @brief Preconditioner ? */
    virtual void
    prec (BlockVector<double> & dst,
          const BlockVector<double> & src) const = 0;

    /** @brief set the Preconditioner */
    virtual void
    set_prec (std::shared_ptr<Equation<dim> > & prec) const = 0;

    virtual void
    solve (BlockVector<double> &dst,
           const BlockVector<double> &src,
           const double tol_default = -1.,
           const unsigned int max_it_default = 0) const = 0;

    /** @brief access to the scattering method. */
    /*virtual void
    scatt_add (BlockVector<double> &dst,
               const BlockVector<double> &src,
               const unsigned int g,
               const unsigned int h) const = 0;*/

    /** @brief access to the up-scattering method. */
    virtual void
    up_add (BlockVector<double> &dst,
            const SnVector &src,
            const unsigned int g) const = 0;

    /** @brief access to the down-scattering method. */
    virtual void
    down_add (BlockVector<double> &dst,
              const SnVector &src,
              const unsigned int g) const = 0;

    /** @brief access to the fission method. */
    /*virtual void
    update_fission_source (SnVector &dst,
                           const SnVector &src) = 0; */

    /** @brief access to the fission method. */
    virtual void
    calculate_fission_density (SnVector &dst,
                            const SnVector &src) = 0;

    /** @brief access to the fission method. */
    virtual void
    distribute_fission_rate (SnVector &dst,
                             const SnVector &src) = 0;

    /** @brief Set the group we are going to perform the group iteration to. */
    virtual void
    set_group (const unsigned int g) = 0;

    /** @brief Is the system symmetric?. */
    virtual bool
    is_symmetric () = 0;

    /** @todo document. */
    virtual void
    set_rhs (BlockVector<double> & dst,
             const BlockVector<double> &src) const = 0;

    /** @todo document. */
    virtual void
    get_angular (SnVector & dst,
                 const SnVector & src,
                 const double keff) const = 0;

    virtual unsigned int
    get_n_angles () const = 0;

    virtual unsigned int
    get_n_groups() const = 0;

    /** @brief Virtual destructor. */
    virtual
    ~Equation ()
    {
    }

    /** @brief Memory consumption. */
    virtual
    double
    memory_consumption () const
    {
      return 0.;
    }

    /** @brief Return the total number of sweeps within energy groups. */
    virtual unsigned int
    get_tot_wg_it () const = 0;

    //virtual bool m_external_preconditioner = 0;

  public:
    virtual bool no_upscatt(const unsigned int g) const = 0;

  };

} // end of namespace Forest

#endif /* FOREST_EQUATION_H */
