/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    algebra/eigen_standard.h
 * @brief   EigenStandard class template declarations
 */

#ifndef FOREST_EIGEN_STANDARD_H
#define FOREST_EIGEN_STANDARD_H

#include "neutronics/state_fwd.h"       // For class State&
#include "algebra/equation_fwd.h"       // for Equation &
#include "algebra/solver_eq.h"
#include "algebra/eigen_prob.h"
#include "algebra/sn_vector.h"          // for SnVector

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <memory>

namespace Forest
{
  using namespace dealii;

  /**
   * @class EigenStandard
   * @ingroup ForestAlgebra
   * @brief Eigenvalue problem in standard form
   *
    @details For the generalized eigenvalue problem

    \f{align}{
     (DLM + S) \phi  = \frac{1}{\lambda} X F^{T} \phi ,
    \f}

    it is rearranged in a standard eigenvalue problem as

    \f{align}{
    A \phi = \lambda \phi ,
    \quad \text{where}
    \quad A :=  (DLM + S)^{-1} X F^{T} .
    \f}

    Here we also store some temporary vectors that we need when performing
    the different operations. This temporary vectors would be for
      - @p m_phi_tmp
      - @p m_fission_density_tmp
      - @p m_fission_source_tmp

   *
   *
   */
  template <int dim>
  class EigenStandard : public EigenProb<dim>
  {
  public:

    /**
     * @brief Constructor of the eigenproblem
     * @param state
     * @param equation
     */
    EigenStandard (State<dim> & state,
                   std::shared_ptr<Equation<dim> > & equation)
      :
        EigenProb<dim>(state, equation)
    {}

    /** @brief Destructor. */
    virtual
    ~EigenStandard ()
    {
    }

    /* Declaration of inherited pure virtual methods */

    /**
    @brief Application of matrix to vector src, and write result into dst.

    @details We perform the matrix vector product with the matrix
    \f{align}{
    A \phi = \lambda \phi ,
    \quad \text{where}
    \quad A :=  (DLM + S)^{-1} X F^{T} .
    \f}

    The value of @p m_phi_tmp has the values from last call to this function
    (or zero for the first call).

    */
    virtual void
    vmult (SnVector &dst,
           const SnVector &src) const;

    /**
     * @brief @copybrief EigenProb::Tvmult()
     * @note Currently NOT implemented for EigenStandard
     */
    virtual void
    Tvmult (SnVector & /* dst */,
           const SnVector & /* src */) const
    {
      Assert(false, ExcMessage("Tvmult not implemented."));
    }

    virtual void
    set_solution_sizes (SnVector & vec) const;

    virtual void
    set_initial_guess () const;

    virtual SnVector &
    get_initial_guess () const;

    /**
     @brief We update the value of the vectors using the solution.

     @details We provide the solution and the eigenvalue, f and keff,
     and we we use the definitions and the relationship by the eigenvalue
     problem to recover intermediate vectors. From

     \f{align}{
     (DLM + S)^{-1} X F^{T} \phi = \lambda \phi ,
     \f}

     we obtain

     \f{align}{
     f := F^{T} \phi
     \f}

     and

     \f{align}{
     q := X f .
     \f}
    */
    virtual void
    update_vectors(const SnVector & vec, const double keff) const;

    virtual double
    memory_consumption () const;

  private:

    /** @brief reference to the State object provided as argument. */
    using EigenProb<dim>::mp_state;
    /** @brief Constant reference to the equation object. */
    using EigenProb<dim>::mp_equation;
    /** @brief Container for the solver. */
    using EigenProb<dim>::m_solver;

  public:
    /** @brief Container for the flux + bv. */
    using EigenProb<dim>::m_phi_tmp;
    /** @brief Container for the fission density. */
    using EigenProb<dim>::m_fission_density_tmp;
    /** @brief Container for the fission source. */
    using EigenProb<dim>::m_fission_source_tmp;
  };

} // end of namespace Forest

#endif /* FOREST_EIGEN_STANDARD_H */
