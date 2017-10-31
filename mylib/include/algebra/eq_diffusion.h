/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eq_diffusion.h
 * @brief   EqDiffusion class template declarations
 */

#ifndef FOREST_EQ_DIFFUSION_H
#define FOREST_EQ_DIFFUSION_H

#include "algebra/equation.h"             // for SnVector, Equation
#include "algebra/solver_eq.h"            // for EqSolverWG and EqPrec
#include "neutronics/state_fwd.h"         // For class State&

#include "deal.II/base/exceptions.h"      // for Assert
#include "deal.II/lac/block_vector.h"     // for BlockVector
#include "deal.II/lac/vector.h"           // for Vector

#include <memory>                         // for shared_ptr
#include <vector>                         // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class ScatteringSource; }
namespace Forest { template <int dim> class FissionSource; }
namespace Forest { template <int dim> class FissionSpectrum; }
namespace Forest { template <int dim> class DiffOp; }
#endif

namespace Forest
{
  using namespace dealii;

  /**
   * @class EqDiffusion
   * @ingroup ForestAlgebra
   * @brief Neutron diffusion equation
   *
   * @details The equation has the following form
   * \f{align*}{
   * \mathcal{C} \Phi = \frac{1}{\lambda} \chi \mathcal{F}\Phi
   * \f}
   * where \f$ \mathcal{C} = \mathcal{L} - \mathcal{S} \f$ is the loss operator
   * composed by the diffusion operator and the scattering term, being the
   * diffusion operator block diagonal matrix defined by
   * \f{align*}{
   * \mathcal{L}_{g,g} \Phi
   * := \nabla \left(\frac{1}{\Sigma_{T}} \nabla \Phi \right) + \Sigma_{T} \Phi
   * \f}
   */
  template <int dim>
  class EqDiffusion : public Equation<dim>
  {
  public:

    /**
     * @brief Constructor for this class
     * @param state
     */
    EqDiffusion (State<dim> & state);

    /** @brief Destructor. */
    virtual
    ~EqDiffusion ()
    {
    }

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_angular_sizes () const;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_flux_bv_sizes () const;

    /** @brief Significant dofs for the incoming flux for every direction. */
    virtual const std::vector<unsigned int> &
    get_fission_density_sizes () const;

    /** @brief Is the system symmetric?. */
    virtual bool
    is_symmetric ()
    {
      return true;
    }

    /** @todo document.  */
    virtual void
    vmult (SnVector & /* dst */ ,
           const SnVector & /* src */ ) const
    {
      Assert(false, ExcMessage("vmult not implemented."));
    }

    /**
     * @brief Application of matrix to vector src, and write result into dst.
     * @details Here we perform, for a given energy group @p m_g, one
     * source iteration. In a compact form it can be expressed as follows
     * \f{align*}{
     * & \hat{x} = (\hat{I} - \hat{D} \hat{L}^{-1} \hat{M}\hat{S}) \hat{b}
     * \f}
     * @param dst
     * @param src
     */
    virtual void
    vmult (BlockVector<double> &dst,
           const BlockVector<double> &src) const;

    /* @brief access to the scattering method. */
    // virtual void
    // scatt_add (BlockVector<double> &dst,
    //            const BlockVector<double> &src,
    //            const unsigned int g,
    //            const unsigned int h) const;

    /** @brief access to the up-scattering method. */
    virtual void
    up_add (BlockVector<double> &dst,
                  const SnVector &src,
                  const unsigned int g) const;

    /** @brief access to the down-scattering method. */
    virtual void
    down_add (BlockVector<double> &dst,
                    const SnVector &src,
                    const unsigned int g) const;

    /* @brief access to the fission method. */
    // virtual void
    // update_fission_source (SnVector &dst,
    //                        const SnVector &src);

    /** @brief access to the fission method. */
    virtual void
    calculate_fission_density (SnVector &dst,
                         const SnVector &src);

    /** @brief access to the fission method. */
    virtual void
    distribute_fission_rate (SnVector &dst,
                             const SnVector &src);

    /**
     * @brief Application of the transpose of a matrix to a vector.
     * Used by some iterative methods.
     */
    virtual void
    Tvmult (SnVector & /* dst */,
            const SnVector & /* src */) const
    {
      Assert(false, ExcNotImplemented());
    }

    /**
     * @brief Application of the transpose of a matrix to a vector.
     * Used by some iterative methods.
     */
    virtual void
    Tvmult (BlockVector<double> & /* dst */,
            const BlockVector<double> & /* src */) const
    {
      Assert(false, ExcNotImplemented());
    }

    /** @brief Preconditioner ? */
    virtual void
    prec (BlockVector<double> & dst,
          const BlockVector<double> & src) const;

    /** @brief set the Preconditioner */
    virtual void
    set_prec (std::shared_ptr<Equation<dim> > & prec) const;

    virtual void
    solve (BlockVector<double> &dst,
           const BlockVector<double> &src,
           const double tol_default = -1.,
           const unsigned int max_it_default = 0) const;

    /**
     * @brief Preparing the right hand side for the source iteration.
     * @details For a given energy group @p m_g, we generate the right
     * hand side as follows
     *
     * \f{align*}{
     * \left[ \begin{array}{c}
     * x_{0} \\
     *  x_{S}
     * \end{array} \right]
     * =
     * \left[ \begin{array}{cc}
     * D & 0 \\
     *  0 & I
     * \end{array} \right]
     * \left[ \begin{array}{c}
     * I \\
     *  P^{T}
     * \end{array} \right]
     * L_{0}^{-1} M Q
     * \f}
     *
     * @param dst
     * @param src
     */
    virtual void
    set_rhs (BlockVector<double> & dst,
             const BlockVector<double> &src) const;


    virtual void
    get_angular (SnVector & dst,
                 const SnVector & src,
                 const double keff) const;

    /** @brief Set the group we are going to perform the group iteration to. */
    virtual void
    set_group (const unsigned int g)
    {
      m_g = g;
    }

    /** @brief Return memory consumption. */
    virtual double
    memory_consumption () const;

    /** @brief Set the group we are going to perform the group iteration to. */
    virtual unsigned int
    get_n_angles () const;

    virtual unsigned int
    get_n_groups() const;

    /** @brief Return the total number of sweeps within energy groups. */
    unsigned int
    get_tot_wg_it () const
    {
      return m_solver.get_tot_wg_it ();
    }

  private:

    /**
     * @brief Application of matrix to vector src, and write result into dst.
     * @details Here we perform, for a given energy group @p m_g, one
     * source iteration. In a matrix form it can be expressed as follows
     *
     * \f{align*}{
     * & \hat{x} = \hat{D} \hat{L}^{-1} \hat{M}\hat{S} \hat{b}
     * \f}
     *
     * We use the following notation for the arguments:
     *
     * \f{align*}
     * \text{src} := b \equiv
     * \left[ \begin{array}{c} b_{0} \\ b_{S} \end{array} \right]
     * \in \mathbb{R}^{I \times M} \times \mathbb{R}^{S} ,
     * \quad \text{and} \quad
     * \text{dst} := x \equiv
     * \left[ \begin{array}{c} x_{0} \\ x_{S} \end{array} \right]
     * \in \mathbb{R}^{I \times M} \times \mathbb{R}^{S}
     * \f}
     *
     * Where \f$ I \f$ is the number of degrees of freedom for the spatial
     * discretization over the mesh, \f$ M \f$ is the number of moments in
     * the scattering operator, and \f$ S \f$ is the number of significant
     * degrees of freedom for the boundary conditions.
     *
     * We perform the following steps:
     *
     * - The first thing is to initialize the vector \f$ x \f$ to zero
     * \f{align*}
     * x :=
     * \left[ \begin{array}{c} x_{0} \\ x_{S} \end{array} \right] =
     * \left[ \begin{array}{c} 0     \\ 0 \end{array} \right]
     * \f}
     *
     * - Then we put the scattering generated into m_rhs
     * \f{align*}
     * v_{1} = S b_{0}
     * \f}
     *
     * - Now, for every direction @p n we do the following operations
     * \f{align*}
     * & \text{for n = 1:N} &  \\
     * &&& v_{2} = M_{n} v_{1} \\
     * &&& v_{2} = v_{2} - P L_{R,n} b_{S,n} \\
     * &&& v_{3} = L_0^{-1} v_{2} \\
     * &&& x_{S} = x_{S} + P^{T} v_{3} \\
     * &&& x_{0} = x_{0} + D_{n} v_{3} \\
     * & \text{endfor} &
     * \f}

     * @param dst
     * @param src
     */
    void
    DiffS (BlockVector<double> &dst,
           const BlockVector<double> &src) const;

    /**
     * @brief parallel version of DL_invMS
     */
    void
    ParDiffS (BlockVector<double> &dst,
              const BlockVector<double> &src) const;


    /** @brief reference to the State object provided as argument. */
    Forest::State<dim> & mp_state;
    /** @brief reference to solver used internally to provide the .solve() method . */
    mutable EqSolverWG<dim> m_solver;
    /** @brief Do we want to do something else than matrix free? */
    bool m_matrix_free;
    /** @brief Transport operator. */
    std::shared_ptr<DiffOp<dim> > m_Diff;
    /** @brief Fission spectrum. */
    std::shared_ptr<FissionSpectrum<dim> > m_Chi;
    /** @brief Fission source. */
    std::shared_ptr<FissionSource<dim> > m_XF;
    /** @brief Scattering source.  */
    std::shared_ptr<ScatteringSource<dim> > m_S;

    /** @brief Preconditioner. */
    mutable std::shared_ptr<Equation<dim> > m_prec;

    mutable unsigned int m_g;
    mutable bool m_external_preconditioner;
    /**
     * @brief Vector containing the Spherical Harmonics moments of the flux.
     * @details The size of this vector is \f$ I \times L \f$, i.e., the number
     * of elements for the spatial mesh times the number of Legendre moments.
     */
    mutable BlockVector<double> m_v1;
    mutable Vector<double> m_v2, m_v3;
    mutable BlockVector<double> m_pv2, m_pv3;

    mutable std::vector<unsigned int> m_angular_sizes;
    mutable std::vector<unsigned int> m_flux_bv_sizes;
    mutable std::vector<unsigned int> m_fission_density_sizes;

  public:
    virtual bool no_upscatt(const unsigned int g) const;
  };

} // end of namespace Forest

#endif /* FOREST_EQ_DIFFUSION_H */
