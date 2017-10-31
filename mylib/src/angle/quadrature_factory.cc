/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   angle/quadrature_factory.cc
 * @brief  Implementation of class template QuadratureFactory
 */

#include "angle/quadrature_factory.h"

#include "angle/quadrature_base.h"
#include "angle/quadrature_gausslegendre.h"
#include "angle/quadrature_chebyshev_legendre.h"
#include "angle/quadrature_levelsymmetric_type1.h"
#include "angle/quadrature_levelsymmetric_type2.h"

#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <iostream>
#include <string>

namespace Forest
{
  //-------------------------------------------------------------------------//
  // 1D QUADRATURES
  //-------------------------------------------------------------------------//
  template <>
  QuadraturePtr<1>
  QuadratureFactory<1>::build (const std::vector<unsigned int> &o,
                               const std::string & type)
  {
    /* The quadrature order for each quadrant is the sn order divided by 2. */
    const unsigned int p_order = o[0]/2;
    /* It should be different from zero with the convention we are using. */
    Assert(p_order != 0, dealii::ExcMessage("polar order equal to zero."));
    /* Generating the Quadrature by calling the derived class constructor. */
    if (type == "GaussLegendre")
      return QuadraturePtr<1>(new QuadratureGaussLegendre<1> (p_order));
    else
      std::cout << "Wrong QuadratureType. GaussLegendre is returned.\n";
    /* This is a default. Think if you want to have many defaults values, or
     * if it is better to give an error. Advertise the default properly if it
     * is used! */
    return QuadraturePtr<1>(new QuadratureGaussLegendre<1> (p_order));
  }

  //-------------------------------------------------------------------------//
  // PURE 2D/3D QUADRATURES
  //-------------------------------------------------------------------------//
  template <int dim>
  QuadraturePtr<dim>
  QuadratureFactory<dim>::build (const std::vector<unsigned int> &o,
                                 const std::string & type)
  {
    if (type == "LevelSymType1")
      return QuadraturePtr<dim>(new QuadratureLevelSymType1<dim> (o[0]));
    else if (type == "LevelSymType2")
      return QuadraturePtr<dim>(new QuadratureLevelSymType2<dim> (o[0]));
    else if (type == "ChebyshevLegendre")
      return QuadraturePtr<dim>(new QuadratureChebyshevLegendre<dim> (o[0], o[1]));
    else
    {
      /** @todo Inconsistent */
      std::cout << "Quadrature type not found. " << std::endl
                << "ChebyshevLegendre will be returned." << std::endl;
    }
    return QuadraturePtr<dim>(new QuadratureLevelSymType2<dim> (o[0]));
  }

  template <int dim>
  void
  QuadratureFactory<dim>::help ()
  {
    std::cout << "The quadratures available in " << dim << "D are:"
        << std::endl;
    switch (dim)
    {
      case 1:
        std::cout << "- gausslegendre" << std::endl;
        break;
      case 2:
        std::cout << "- levelsymmetric_type1" << std::endl
        << "- levelsymmetric_type2" << std::endl
        << "- ChebyshevLegendre" << std::endl;
        break;
      case 3:
        std::cout << "- levelsymmetric_type1" << std::endl
          << "- levelsymmetric_type2" << std::endl
          << "- ChebyshevLegendre" << std::endl;
        break;
      default:
        std::cout << "Dimension not allow";
        break;
    }
    std::cout << std::endl;
  }

  // Explicit instantiations
  template class QuadratureFactory<1> ;
  template class QuadratureFactory<2> ;
  template class QuadratureFactory<3> ;

} // end of namespace Forest
