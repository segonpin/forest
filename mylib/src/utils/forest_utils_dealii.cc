/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   utils/forest_utils_dealii.cc
 * @brief  Implementation of some functions
 */

#include "utils/forest_utils_dealii.h"

#include <iostream>

namespace Forest
{
  //-------------------------------------------------------------------------
  // template function definition

  /**
   * @brief 1D specialization when the shape is not provided
   * (so it is assumed that the reactor fills the whole box)
   */
  template <typename Number>
  dealii::Table<1, Number>
  shape_to_table_1 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials)
  {
    // Defining the table and filling with (-1)
    dealii::Table<1, Number> existing_cells (n_nodes[0]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int i = 0; i < n_nodes[0]; ++i)
    {
      existing_cells[i] = materials[mat];
      ++mat;
    }
    return existing_cells;
  }

  // template function definition
  /**
   * @brief 2D specialization when the shape is not provided
   * (so it is assumed that the reactor fills the whole box)
   */
  template <typename Number>
  dealii::Table<2, Number>
  shape_to_table_2 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials)
  {
    // Defining the table and filling with (-1)
    dealii::Table<2, Number> existing_cells (n_nodes[0], n_nodes[1]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int j = 0; j < n_nodes[1]; ++j)
    {
      for (unsigned int i = 0; i < n_nodes[0]; ++i)
      {
        existing_cells[i][j] = materials[mat];
        ++mat;
      }
    }
    return existing_cells;
  }

  // template function definition
  /**
   * @brief 3D specialization when the shape is not provided
   * (so it is assumed that the reactor fills the whole box)
   */
  template <typename Number>
  dealii::Table<3, Number>
  shape_to_table_3 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials)
  {
    // Defining the table and filling with (-1)
    dealii::Table<3, Number> existing_cells (n_nodes[0], n_nodes[1],
        n_nodes[2]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int k = 0; k < n_nodes[2]; ++k)
    {
      for (unsigned int j = 0; j < n_nodes[1]; ++j)
      {
        for (unsigned int i = 0; i < n_nodes[0]; ++i)
        {
          existing_cells[i][j][k] = materials[mat];
          ++mat;
        }
      }
    }
    return existing_cells;
  }

  //-------------------------------------------------------------------------
  // Specialization for dim = 1,  Number = unsigned char
  /**
   * @brief specialization
   */
  template <>
  dealii::Table<1, unsigned char>
  shape_to_table<1, unsigned char> (std::vector<unsigned int> &n_nodes,
                                    std::vector<unsigned char> &materials)
  {
    return shape_to_table_1<unsigned char> (n_nodes, materials);
  }

  // Specialization for dim = 2,  Number = unsigned char
  /**
   * @brief specialization
   */
  template <>
  dealii::Table<2, unsigned char>
  shape_to_table<2, unsigned char> (std::vector<unsigned int> &n_nodes,
                                    std::vector<unsigned char> &materials)
  {
    return shape_to_table_2<unsigned char> (n_nodes, materials);
  }

  // Specialization for dim = 3,  Number = unsigned char
  /**
   * @brief specialization
   */
  template <>
  dealii::Table<3, unsigned char>
  shape_to_table<3, unsigned char> (std::vector<unsigned int> &n_nodes,
                                    std::vector<unsigned char> &materials)
  {
    return shape_to_table_3<unsigned char> (n_nodes, materials);
  }

  //-------------------------------------------------------------------------

  /**
   * @brief 1D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  dealii::Table<1, Number>
  shape_to_table_1 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> & /* geom_shape */)
  {
    // Defining the table and filling with (-1)
    dealii::Table<1, Number> existing_cells (n_nodes[0]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int i = 0; i < n_nodes[0]; ++i)
    {
      existing_cells[i] = materials[mat];
      ++mat;
    }
    return existing_cells;
  }

  /**
   * @brief 2D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  dealii::Table<2, Number>
  shape_to_table_2 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> &geom_shape)
  {
    // Defining the table and filling with (-1)
    dealii::Table<2, Number> existing_cells (n_nodes[0], n_nodes[1]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int j = 0; j < n_nodes[1]; ++j)
    {
      for (unsigned int i = geom_shape[0].x_begin[j];
          i < geom_shape[0].x_end[j]; ++i)
      {
        existing_cells[i][j] = materials[mat];
        ++mat;
      }
    }
    return existing_cells;
  }

  /**
   * @brief 3D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  dealii::Table<3, Number>
  shape_to_table_3 (std::vector<unsigned int> &n_nodes,
                    std::vector<Number> &materials,
                    std::vector<xy_shape> &geom_shape)
  {
    // Defining the table and filling with (-1)
    dealii::Table<3, Number> existing_cells (n_nodes[0], n_nodes[1],
        n_nodes[2]);
    existing_cells.fill ((Number) (-1));

    // Filling the table with the values
    unsigned int mat = 0;
    for (unsigned int k = 0; k < n_nodes[2]; ++k)
    {
      for (unsigned int j = 0; j < n_nodes[1]; ++j)
      {
        for (unsigned int i = geom_shape[k].x_begin[j];
            i < geom_shape[k].x_end[j]; ++i)
        {
          existing_cells[i][j][k] = materials[mat];
          ++mat;
        }
      }
    }
    return existing_cells;
  }

  //-------------------------------------------------------------------------

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  dealii::Table<1, unsigned char>
  shape_to_table (std::vector<unsigned int> &n_nodes,
                  std::vector<unsigned char> &materials,
                  std::vector<xy_shape> &flat_shape)
  {
    return shape_to_table_1 (n_nodes, materials, flat_shape);
  }

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  dealii::Table<2, unsigned char>
  shape_to_table (std::vector<unsigned int> &n_nodes,
                  std::vector<unsigned char> &materials,
                  std::vector<xy_shape> &flat_shape)
  {
    return shape_to_table_2 (n_nodes, materials, flat_shape);
  }

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  dealii::Table<3, unsigned char>
  shape_to_table (std::vector<unsigned int> &n_nodes,
                  std::vector<unsigned char> &materials,
                  std::vector<xy_shape> &flat_shape)
  {
    return shape_to_table_3 (n_nodes, materials, flat_shape);
  }

  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------

  /**
   * @brief 1D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  void
  print_table_1 (dealii::Table<1, Number> &table)
  {
    std::cout << "printing table \n";
    for (unsigned int i = 0; i < table.size (0); ++i)
    {
      std::cout << static_cast<double> (table[i]) << " ";
    }
    std::cout << std::endl;
  }

  /**
   * @brief 2D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  void
  print_table_2 (dealii::Table<2, Number> &table)
  {
    std::cout << "printing table \n";
    for (unsigned int j = 0; j < table.size (1); ++j)
    {
      for (unsigned int i = 0; i < table.size (0); ++i)
      {
        std::cout << static_cast<double> (table[i][j]) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  /**
   * @brief 3D specialization when the shape is provided
   * (so it is assumed that the reactor does not fill the whole box)
   */
  template <typename Number>
  void
  print_table_3 (dealii::Table<3, Number> &table)
  {
    std::cout << "printing table \n";
    for (unsigned int k = 0; k < table.size (2); ++k)
    {
      for (unsigned int j = 0; j < table.size (1); ++j)
      {
        for (unsigned int i = 0; i < table.size (0); ++i)
        {
          std::cout << static_cast<double> (table[i][j][k]) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
    }
    std::cout << std::endl;
  }

  //-------------------------------------------------------------------------

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  void
  print_table (dealii::Table<1, unsigned char> &table)
  {
    return print_table_1 (table);
  }

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  void
  print_table (dealii::Table<2, unsigned char> &table)
  {
    return print_table_2 (table);
  }

  /**
   * @brief This should NOT be used, only defined to specialization below:
   */
  template <>
  void
  print_table (dealii::Table<3, unsigned char> &table)
  {
    return print_table_3 (table);
  }

} // end of namespace Forest
