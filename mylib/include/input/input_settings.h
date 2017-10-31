/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input_settings.h
 * @brief  InputSettings class template declarations
 */

#ifndef FOREST_INPUT_SETTINGS_H
#define FOREST_INPUT_SETTINGS_H

#include <string>
#include <vector>

// This is for the xml parser from boost
#include <boost/property_tree/ptree_fwd.hpp>

namespace Forest
{
  /**
   @brief Class containing the input and output files.
   @details Information obtained/generated mainly from the settings file.
   @ingroup ForestInput
   */
  class InputFiles
  {
  public:

    /**
     @brief Constructor initializing all the members to empty strings.
     */
    InputFiles ();

  protected:

    std::string geom; /**< Input geometry file address. */
    std::string mat; /**< Input materials file address. */
    std::string out; /**< Output file address. */
    std::string mesheps; /**< Output mesh and Geometry file address. */
    std::string meshvtk; /**< Output mesh and Geometry file address. */
    std::string vtk; /**< Output vtk file address. */
    std::string bcs; /**< Output/input Boundary conditions file address. */
    std::string log; /**< Output-log file address. */
    bool bcs_from_file; /**< Boundary conditions  from external file? */

    bool m_initialized;

    void
    set (const std::string &path,
         const std::string &filename_no_ext);

    /**
     @brief Read the names of the input files from a xml tree, and save it.
     @param pt
     @param path
     @param filename_no_ext
     */
    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);

    /**
     @brief Write the content of the files in a xml tree.
     @param pt
     */
    void
    write (boost::property_tree::ptree & pt) const;

    /** To have access to protected members. */
    friend class InputSettings;

  };

  /**
   @brief Class containing input information for the coarsening
   @details Information obtained/generated mainly from the settings file.
   */
  /*class InputCoarseningBcs
  {
  public:
    InputCoarseningBcs ();
  protected:
    unsigned int spatial;
    unsigned int angular;
    unsigned int energy;
    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);
    void
    write (boost::property_tree::ptree & pt) const;
    friend class InputSettings;
  };*/

  /**
   @brief Class containing input information for the problem
   @details Information obtained/generated mainly from the settings file.
   */
  class InputProblem
  {
  public:
    InputProblem ();
  protected:
    std::string m_type;
    bool m_use_transport;
    bool m_use_dsa;
    std::string m_quad_name;
    bool m_product_quadrature;
    std::vector<unsigned int> m_order;

    void
    set (const std::string type = std::string ("transport"),
         bool product_quadrature = false,
         std::string quad = "LevelSymType2",
         std::vector<unsigned int> order =
         { 2 })
    {
      set_type (type);
      set_quad (product_quadrature, quad, order);
    }
    void
    set_type (const std::string type);
    void
    set_quad (bool product_quadrature = false,
              std::string quad = "LevelSymType2",
              std::vector<unsigned int> order =
              { 2 });

    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);
    void
    write (boost::property_tree::ptree & pt) const;

    /** To have access to protected members. */
    friend class InputSettings;
  };

  /**
   @brief Class containing input information for the algebra
   @details Information obtained/generated mainly from the settings file.
   */
  class InputAlgebra
  {
  public:
    InputAlgebra ();
  protected:
    bool m_matrix_free;
    bool m_use_power_iteration;
    unsigned int m_n_eigenvalues;
    unsigned int m_n_cv;
    bool m_use_fission_density;

    std::string m_eig_solver;
    double m_eig_tol;
    unsigned int m_eig_max_it;

    std::string m_mg_solver;
    double m_mg_tol;
    unsigned int m_mg_max_it;

    std::string m_inner_solver;
    double m_inner_tol;
    unsigned int m_inner_max_it;
    void
    set (bool matrix_free = true,
         const bool use_fission_density = true,
         const std::string eig_solver = "PI",
         const double eig_tol = 1.e-6,
         const unsigned int eig_max_it = 1000,
         const std::string mg_solver = "GS",
         const double mg_tol = 1.e-6,
         const unsigned int mg_max_it = 1000,
         const std::string inner_solver = "Krylov",
         const double inner_tol = 1.e-7,
         const unsigned int inner_max_it = 1000);

    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);
    void
    write (boost::property_tree::ptree & pt) const;

    /** To have access to protected members. */
    friend class InputSettings;
  };

  /**
   @brief Class containing input information for the Finite elements
   @details Information obtained/generated mainly from the settings file.
   */
  class InputFESettings
  {
  public:
    InputFESettings ();
  protected:
    unsigned int m_degree; /**< Degree of the finite element expansion. */
    unsigned int m_n_ref; /**< Number of global refinements to be applied. */

    void
    set (unsigned int degree = 0,
         unsigned int n_ref = 0);

    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);
    void
    write (boost::property_tree::ptree & pt) const;

    /** To have access to protected members. */
    friend class InputSettings;
  };

  /**
   @brief Class containing output information for the problem
   */
  class InputOutput
  {
  public:
    InputOutput ();
  protected:
    bool homogenize;
    bool print_mesh;
    bool vtk_power;
    bool vtk_flux;
    bool vtk_mat;
    bool ass_bcs;

    void
    set (const std::vector<bool> & options);

    void
    read_xml (const boost::property_tree::ptree & pt,
              const std::string &path,
              const std::string &filename_no_ext);
    void
    write (boost::property_tree::ptree & pt) const;

    /** To have access to protected members. */
    friend class InputSettings;
  };

  /**
   * @brief Container for the settings of the problem.
   */
  class InputSettings
  {

  public:

    /* @brief Empty constructor */
    InputSettings ();

    /**
     @brief Initialize the classes member and read the settings.xml file.
     @details All the work is done inside the load and save functions before.
     @param path_plus_file_plus_ext_
     */
    InputSettings (const std::string & path_plus_file_plus_ext_);

    /**
     * @brief Once the property tree is created, extract all the information.
     * @param path
     * @param filename_no_ext
     */
    void
    load_xml (const std::string &path,
              const std::string &filename_no_ext);

    /**
     * @brief This performs the inverse operation of the @ref load function,
     *        printing the information in a file.
     * @param path
     * @param filename_no_ext
     */
    void
    save_xml (const std::string &path,
              const std::string &filename_no_ext);

  public:

    std::string settings_ext = ".settings.xml"; /**< Extension used. */
    std::string path_plus_file_plus_ext; /**< Full path input file. */
    std::string path; /**< Path for the problem. */
    std::string problem_name; /**< Problem name provided. */

    InputFiles files; /**< Settings for the files. */
    //InputCoarseningBcs coarsening_bcs; /**< Settings for boundary conditions. */
    InputProblem problem; /**< Settings for the problem to be solved. */
    InputAlgebra algebra; /**< Settings for the algebra. */
    InputFESettings fe_settings; /**< Settings for finite elements. */
    InputOutput output; /**< Settings for the output. */

  public:

    void
    set_files (std::string path,
               std::string problem_name)
    {
      files.set (path, problem_name);
    }

    void
    set_problem (std::string path,
                 std::string problem_name)
    {
      files.set (path, problem_name);
    }

    std::string
    get_path_plus_file () const
    {
      return path + problem_name;
    }

    const std::string &
    get_path_plus_file_plus_ext () const
    {
      return path_plus_file_plus_ext;
    }

    const std::string &
    get_file_geom () const
    {
      return files.geom;
    }

    const std::string &
    get_file_mat () const
    {
      return files.mat;
    }

    const std::string &
    get_file_out () const
    {
      return files.out;
    }

    const std::string &
    get_file_mesheps () const
    {
      return files.mesheps;
    }

    const std::string &
    get_file_meshvtk () const
    {
      return files.meshvtk;
    }

    const std::string &
    get_file_vtk () const
    {
      return files.vtk;
    }

    const std::string &
    get_file_bcs () const
    {
      return files.bcs;
    }

    std::string
    get_log () const
    {
      return path + files.log;
    }

    /**
     *
     * @return
     */
    /*unsigned int
    get_coarse_spatial () const
    {
      return coarsening_bcs.spatial;
    }*/

    /**
     *
     * @return
     */
    /*unsigned int
    get_coarse_angular () const
    {
      return coarsening_bcs.angular;
    }*/

    /**
     *
     * @return
     */
    /*unsigned int
    get_coarse_energy () const
    {
      return coarsening_bcs.energy;
    }*/

    bool
    use_transport () const
    {
      return problem.m_use_transport;
    }

    bool
    use_dsa () const
    {
      return problem.m_use_dsa;
    }

    std::vector<unsigned int>
    get_sn_order () const
    {
      return problem.m_order;
    }

    std::string
    get_quad_name () const
    {
      return problem.m_quad_name;
    }

    bool
    get_matrix_free () const
    {
      return algebra.m_matrix_free;
    }

    bool
    use_power_iteration () const
    {
      return algebra.m_use_power_iteration;
    }

    unsigned int
    get_n_eigenvalues () const
    {
      return algebra.m_n_eigenvalues;
    }

    unsigned int
    get_n_column_vectors () const
    {
      return algebra.m_n_cv;
    }

    bool
    use_fission_density () const
    {
      return algebra.m_use_fission_density;
    }

    std::string
    get_eig_solver () const
    {
      return algebra.m_eig_solver;
    }

    double
    get_eig_tol () const
    {
      return algebra.m_eig_tol;
    }

    unsigned int
    get_eig_max_it () const
    {
      return algebra.m_eig_max_it;
    }

    std::string
    get_mg_solver () const
    {
      return algebra.m_mg_solver;
    }

    double
    get_mg_tol () const
    {
      return algebra.m_mg_tol;
    }

    unsigned int
    get_mg_max_it () const
    {
      return algebra.m_mg_max_it;
    }

    std::string
    get_inner_solver () const
    {
      return algebra.m_inner_solver;
    }

    double
    get_inner_tol () const
    {
      return algebra.m_inner_tol;
    }

    unsigned int
    get_inner_max_it () const
    {
      return algebra.m_inner_max_it;
    }

    unsigned int
    get_n_ref () const
    {
      return fe_settings.m_n_ref;
    }

    unsigned int
    get_fe_degree () const
    {
      return fe_settings.m_degree;
    }

    bool
    get_homogenize_flag () const
    {
      return output.homogenize;
    }

    bool
    get_print_mesh_flag () const
    {
      return output.print_mesh;
    }

    bool
    get_print_vtk_flag () const
    {
      return (output.vtk_mat or output.vtk_power or output.vtk_flux);
    }

    bool
    get_print_vtk_mat_flag () const
    {
      return output.vtk_mat;
    }

    bool
    get_print_vtk_power_flag () const
    {
      return output.vtk_power;
    }

    bool
    get_print_vtk_flux_flag () const
    {
      return output.vtk_flux;
    }

    bool
    get_print_ass_bcs_flag () const
    {
      return output.ass_bcs;
    }

    bool
    get_bcs_from_file () const
    {
      return files.bcs_from_file;
    }

    void
    set_input_files (const std::string &path,
                     const std::string &filename_no_ext)
    {
      files.set (path, filename_no_ext);
    }

    void
    set_problem (const std::string type = std::string ("transport"),
                 bool product_quadrature = false,
                 std::string quad = "LevelSymType2",
                 std::vector<unsigned int> order =
                 { 2 })
    {
      problem.set (type, product_quadrature, quad, order);
    }

    void
    set_algebra (bool matrix_free,
                 const bool use_fission_density,
                 const std::string eig_solver,
                 const double eig_tol,
                 const unsigned int eig_max_it,
                 const std::string mg_solver,
                 const double mg_tol,
                 const unsigned int mg_max_it,
                 const std::string inner_solver,
                 const double inner_tol,
                 const unsigned int inner_max_it)
    {
      algebra.set (matrix_free, use_fission_density,
          eig_solver, eig_tol, eig_max_it,
          mg_solver, mg_tol, mg_max_it,
          inner_solver, inner_tol, inner_max_it);
    }

    void
    set_fe_settings (unsigned int degree,
                     unsigned int n_ref)
    {
      fe_settings.set (degree, n_ref);
    }

    void
    set_output (const std::vector<bool> & options)
    {
      output.set (options);
    }

  };

} // end of namespace Forest
#endif /* FOREST_INPUT_SETTINGS_H */
