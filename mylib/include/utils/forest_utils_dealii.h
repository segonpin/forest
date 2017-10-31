#ifndef FOREST_UTILS_DEALII_H
#define FOREST_UTILS_DEALII_H

#include "utils/forest_utils_base.h"                        //for xy_shape

#include "deal.II/base/point.h"                             // for Point
#include "deal.II/base/table.h"                             // for Table
#include "deal.II/base/types.h"
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <vector>                                           // for vector

//
// Here we have different functions that can be useful for different small
// tasks. The namespace is documented in forest_utils.h
//

namespace Forest
{
  using namespace dealii;

//-------------------------------------------------------------------------
// File utils

//-------------------------------------------------------------------------
  /**
   * @brief Fill dealii Table<dim, number> with std::vector<std::vector<number> >
   * @todo generalize for more dimensions
   */
  template <int dim>
  Table<dim, types::material_id>
  stdvector_to_table (std::vector<unsigned int> n_nodes,
                      std::vector<std::vector<types::material_id> > materials);

#ifndef DOXYGEN
  template <>
  Table<1, types::material_id>
  stdvector_to_table<1> (
      std::vector<unsigned int> n_nodes,
      std::vector<std::vector<types::material_id> > materials);
#endif
  //-------------------------------------------------------------------------

//   typedef struct
//   {
//     std::vector<unsigned int> x_begin, x_end;
//   } xy_shape;

//-------------------------------------------------------------------------

  template <typename Number>
  Table<1, Number>
  shape_to_table_1 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials);

  template <typename Number>
  Table<2, Number>
  shape_to_table_2 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials);

  template <typename Number>
  Table<3, Number>
  shape_to_table_3 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials);

  //-------------------------------------------------------------------------

  /**
   * @brief fill dealii Table<dim, number> with reactor filling whole box
   */
  template <int dim, typename Number>
  Table<dim, Number>
  shape_to_table (std::vector<unsigned int> &n_nodes,
                  std::vector<Number> &materials);

  //-------------------------------------------------------------------------

  template <typename Number>
  Table<1, Number>
  shape_to_table_1 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> &flat_shape);

  template <typename Number>
  Table<2, Number>
  shape_to_table_2 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> &flat_shape);

  template <typename Number>
  Table<3, Number>
  shape_to_table_3 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> &flat_shape);

  //-------------------------------------------------------------------------

  /**
   * @brief fill dealii Table<dim, number> with reactor filling xy_shape
   */
  template <int dim, typename Number>
  Table<dim, Number>
  shape_to_table (std::vector<unsigned int> &n_nodes,
                  std::vector<Number> &materials,
                  std::vector<xy_shape> &flat_shape);

//-------------------------------------------------------------------------

  /**
   * @brief Print content of dealii Table<1, number>
   */
  template <typename Number>
  void
  print_table_1 (Table<1, Number> &table);

  /**
   * @brief Print content of dealii Table<1, number>
   */
  template <typename Number>
  void
  print_table_2 (Table<2, Number> &table);

  /**
   * @brief Print content of dealii Table<1, number>
   */
  template <typename Number>
  void
  print_table_3 (Table<3, Number> &table);

//-------------------------------------------------------------------------

  /**
   * @brief Print content of dealii Table<dim, number>
   */
  template <int dim, typename Number>
  void
  print_table (Table<dim, Number> &table);

//-------------------------------------------------------------------------

  /**
   * @todo Document vector_to_point
   */
  template <int dim>
  inline Point<dim>
  vector_to_point (const std::vector<double> & v)
  {
    const Point<dim> beta =
        dim == 1 ? Point<dim> (v[0]) : dim == 2 ? Point<dim> (v[0], v[1]) :
        dim == 3 ? Point<dim> (v[0], v[1], v[2]) : Point<dim> ();
    return beta;
  }

  /**
   * @brief We create a new Point<dim> object.
   * @details Depending on the template parameter @p dim we will use the
   * corresponding number of input elements. In case the number of provided
   * elements is small than the dimension the default will be to fill with a
   * zero value.
   * @param v0
   * @param v1
   * @param v2
   * @return
   */
  template <int dim>
  inline Point<dim>
  new_point (const double v0 = 0.0, const double v1 = 0.0,
             const double v2 = 0.0)
  {
    const Point<dim> beta = dim == 1 ? Point<dim> (v0) :
                            dim == 2 ? Point<dim> (v0, v1) :
                            dim == 3 ? Point<dim> (v0, v1, v2) : Point<dim> ();
    return beta;
  }

  /**
   * @details The notation is following the following convention (for dim == 3):
   * @verbatim
   * p_normal(0) = {-1,  0,  0}
   * p_normal(1) = { 1,  0,  0}
   * p_normal(2) = { 0, -1,  0}
   * p_normal(3) = { 0,  1,  0}
   * p_normal(4) = { 0,  0, -1}
   * p_normal(5) = { 0,  0,  1}
   * @endverbatim
   * @param f
   * @return
   */
  template <int dim>
  inline Tensor<1, dim>
  p_normal (unsigned int f)
  {
    return vector_to_point<dim> (r_normal<dim> (f));
  }

  /**
   * @brief Find the face for which the vector is an outgoing vector
   * @details considering the numbering given by p_normal
   * @verbatim
   * p_normal(0) = {-1,  0,  0}
   * p_normal(1) = { 1,  0,  0}
   * p_normal(2) = { 0, -1,  0}
   * p_normal(3) = { 0,  1,  0}
   * p_normal(4) = { 0,  0, -1}
   * p_normal(5) = { 0,  0,  1}
   * @endverbatim
   * we return the index f such that p_normal(f) == normal
   * @param normal
   * @return face
   * @note This is probably very inefficient, but for the moment does not
   * matter because it is done in a post processing. Can be made faster if all
   * the constants already exist, and the function is then a member of the
   * class that uses it.
   */
  template <int dim>
  unsigned int
  normal2face (const Tensor<1, dim> & normal)
  {
    /* We need this constant to define the range */
    const unsigned int faces_per_cell = 2*dim;
    /* Initialization of the face to be returned with an invalid value */
    unsigned int face = faces_per_cell;
    /* Tolerance to decide if the vectors are the same */
    const double tol = 1.e-10;
    /* Run over the faces and compare */
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      /* we check if the vectors are the same direction and orientation */
      if (std::abs (p_normal<dim> (f) * normal - 1.) < tol)
      {
        /* If they are aligned, we store the value and stop the for loop */
        face = f;
        continue;
      }
    }
    /* Verify that the value is inside the expected range */
    Assert(face < faces_per_cell, ExcMessage("face out of bounds"));
    /* Return the value */
    return face;
  }

  template <int dim>
  unsigned int
  normal2face (const Point<dim> & normal)
  {
    /* We need this constant to define the range */
    const unsigned int faces_per_cell = 2*dim;
    /* Initialization of the face to be returned with an invalid value */
    unsigned int face = faces_per_cell;
    /* Tolerance to decide if the vectors are the same */
    const double tol = 1.e-10;
    /* Run over the faces and compare */
    for (unsigned int f = 0; f < faces_per_cell; ++f)
    {
      /* we check if the vectors are the same direction and orientation */
      if (std::abs (p_normal<dim> (f) * normal - 1.) < tol)
      {
        /* If they are aligned, we store the value and stop the for loop */
        face = f;
        continue;
      }
    }
    /* Verify that the value is inside the expected range */
    Assert(face < faces_per_cell, ExcMessage("face out of bounds"));
    /* Return the value */
    return face;
  }


} // end of namespace Forest

#endif
