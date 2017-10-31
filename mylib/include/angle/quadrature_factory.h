/**
 * @author  Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file    angle/quadrature_factory.h
 * @brief   QuadratureFactory class template declarations
 */

#ifndef FOREST_QUADRATURE_FACTORY_H
#define FOREST_QUADRATURE_FACTORY_H

#include "angle/quadrature_base_fwd.h"

#include <memory>
#include <string>
#include <vector>

//
// Classes used for quadrature. The group is documented in forest_angle.h
//

namespace Forest
{
  /** @brief We create an alias to the shared_ptr with template aliasing. */
  template <int dim>
  using QuadraturePtr = std::shared_ptr<QuadratureBase<dim> >;

  /**
   * @class QuadratureFactory
   * @brief Generate the quadratures and return the object.
   * @ingroup ForestAngle
   *
   * The quadrature factory is a class intended to host the functions
   * used to create the quadratures and give them back to a
   * base class object. Thus, this class actually needs no object,
   * and because of this the member functions are declared as static.
   *
   * @note Lets hope the compiler uses the optimization from the 
   * c++11 rvalue ...
   *
   * @todo I have to change the input arguments for the input file,
   * in such a way that I look into the input file to discern which
   * quadrature should be implemented.
   *
   * @todo I have to add more quadrature types.
   */
  template <int dim>
  class QuadratureFactory
  {
  public:

    /**
     * @brief Build the quadrature.
     *
     * A quadrature is generated and the object is returned. The choice
     * of the quadrature is obtained from @p name_, while the order of
     * the quadrature is defined by @p order_. The quadrature is created
     * by a derived class, and the object returned is assigned to the base
     * class (which stores only the data and some members function to
     * access it.
     *
     * @param order_
     * @param type
     * @return
     */

    static QuadraturePtr<dim>
    build (const std::vector<unsigned int> &order_,
           const std::string & type = "ChebyshevLegendre");

    /** @brief Returns a list of available quadratures. */
    static void
    help ();

  };

} // end of namespace Forest

#endif
