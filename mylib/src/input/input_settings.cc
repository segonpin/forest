/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input_settings.cc
 * @brief  Implementation of class InputSettings
 */

#include "input/input_settings.h"

#include "utils/forest_utils_base.h"

// This is for the xml parser from boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION < 105500)
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#endif

namespace Forest
{

  // Input Files ---------------------------------------------------------------

  InputFiles::InputFiles ()
      : geom (""),
        mat (""),
        out (""),
        mesheps (""),
        meshvtk (""),
        vtk (""),
        bcs (""),
        log (""),
        bcs_from_file (false),
        m_initialized(false)
  {

  }

  void
  InputFiles::set (const std::string &path,
                    const std::string &filename_no_ext)
  {
    m_initialized = true;

    geom = path + filename_no_ext + ".geom.xml";
    mat = path + filename_no_ext + ".mat.xml";
    out = path + filename_no_ext + ".out.xml";
    mesheps = path + filename_no_ext + ".mesh.eps";
    meshvtk = path + filename_no_ext + ".mesh.vtk";
    vtk = path + filename_no_ext + ".vtk";
  }

  void
  InputFiles::read_xml (const boost::property_tree::ptree & pt,
                    const std::string &path,
                    const std::string &filename_no_ext)
  {
    m_initialized = true;
    geom = path
        + pt.get<std::string> ("settings.input_files.geom",
            filename_no_ext + ".geom.xml");
    mat = path
        + pt.get<std::string> ("settings.input_files.mat",
            filename_no_ext + ".mat.xml");
    out = path
        + pt.get<std::string> ("settings.input_files.out",
            filename_no_ext + ".out.xml");
    mesheps = path
        + pt.get<std::string> ("settings.input_files.mesheps",
            filename_no_ext + ".mesh.eps");
    meshvtk = path
        + pt.get<std::string> ("settings.input_files.meshvtk",
            filename_no_ext + ".mesh.vtk");
    vtk = path
        + pt.get<std::string> ("settings.input_files.vtk",
            filename_no_ext + ".vtk");

    bcs = pt.get<std::string> ("settings.input_files.bcs", "");
    if (bcs.size () != 0)
      bcs_from_file = true;
    if (bcs_from_file)
      bcs = path + bcs;
  }

  void
  InputFiles::write (boost::property_tree::ptree & pt) const
  {
    pt.add ("settings.<xmlcomment>", "input_files");
    pt.add ("settings.input_files", "");
    pt.put ("settings.input_files.geom", geom);
    pt.put ("settings.input_files.mat", mat);
    pt.put ("settings.input_files.out", out);
    pt.put ("settings.input_files.mesheps", mesheps);
    pt.put ("settings.input_files.meshvtk", meshvtk);
    pt.put ("settings.input_files.vtk", vtk);
    if (bcs_from_file)
      pt.put ("settings.input_files.bcs", bcs);
  }

  // InputCoarseningBcs --------------------------------------------------------
  /*
  InputCoarseningBcs::InputCoarseningBcs ()
      : spatial (0),
        angular (0),
        energy (0)
  {
  }

  void
  InputCoarseningBcs::read_xml (const boost::property_tree::ptree & pt,
                            const std::string & ,
                            const std::string & )
  {
    spatial = pt.get<unsigned int> ("settings.coarsening_bcs.spatial", 0);
    angular = pt.get<unsigned int> ("settings.coarsening_bcs.angular", 0);
    energy = pt.get<unsigned int> ("settings.coarsening_bcs.energy", 0);
  }

  void
  InputCoarseningBcs::write (boost::property_tree::ptree & pt) const
  {
    pt.add ("settings.<xmlcomment>", "coarsening_bcs");
    pt.add ("settings.coarsening_bcs", "");
    pt.put ("settings.coarsening_bcs.spatial", spatial);
    pt.put ("settings.coarsening_bcs.angular", angular);
    pt.put ("settings.coarsening_bcs.energy", energy);
  }*/

  // InputProblem --------------------------------------------------------------

  InputProblem::InputProblem ()
      : m_type ("transport"),
        m_use_transport (true),
        m_use_dsa (false),
        m_quad_name ("LevelSymType2"),
        m_product_quadrature (false),
        m_order (2)
  {
  }

  void
  InputProblem::set_type (std::string type)
  {
    const std::string string_transport ("transport");
    const std::string string_diffusion ("diffusion");
    m_type = type;

    if (m_type.compare (string_transport) == 0)
    {
      m_use_transport = true;
    }
    else if (m_type.compare (string_diffusion) == 0)
    {
      m_use_transport = false;
    }
    else
    {
      m_type = string_transport; // Default or exception?
    }
  }

  void
  InputProblem::set_quad (bool product_quadrature,
      std::string quad,
      std::vector<unsigned int> order)
  {
    m_product_quadrature = product_quadrature;
    m_quad_name = quad;
    m_order = order;
  }

  void
  InputProblem::read_xml (const boost::property_tree::ptree & pt,
                      const std::string & /* path */,
                      const std::string & /* filename_no_ext */)
  {
    m_type = pt.get<std::string> ("settings.problem.<xmlattr>.type");

    const std::string string_transport ("transport");
    const std::string string_diffusion ("diffusion");
    if (m_type.compare (string_transport) == 0)
    {
      m_use_transport = true;
      m_use_dsa = pt.get<bool> ("settings.problem.<xmlattr>.use_dsa", false);

      std::string bin = pt.get<std::string> (
          "settings.problem.quadrature.<xmlattr>.type", "");
      if (bin.compare (std::string ("product")) == 0)
        m_product_quadrature = true;

      if (m_product_quadrature)
      {
        m_quad_name = pt.get<std::string> (
            "settings.problem.quadrature.<xmlattr>.name", "ChebyshevLegendre");
        m_order.resize (2);
        m_order[0] = pt.get<unsigned int> (
            "settings.problem.quadrature.<xmlattr>.azimuthal_order");
        m_order[1] = pt.get<unsigned int> (
            "settings.problem.quadrature.<xmlattr>.polar_order");
      }
      else
      {
        m_quad_name = pt.get<std::string> (
            "settings.problem.quadrature.<xmlattr>.name", "LevelSymType2");
        m_order.resize (1);
        m_order[0] = pt.get<unsigned int> (
            "settings.problem.quadrature.<xmlattr>.sn");
      }
    }
    else if (m_type.compare (string_diffusion) == 0)
    {
      m_use_transport = false;
    }
    else
    {
      //The default is also diffusion
      m_use_transport = false;
    }
  }

  void
  InputProblem::write (boost::property_tree::ptree & pt) const
  {
    pt.add ("settings.<xmlcomment>", "problem");
    pt.add ("settings.problem", "");
    pt.put ("settings.problem.<xmlattr>.type", m_type);
    if (m_use_transport)
    {
      if (m_use_dsa)
      {
        pt.put ("settings.problem.<xmlattr>.use_dsa", m_use_dsa);
      }

      pt.put ("settings.problem.quadrature.name", m_quad_name);
      if (m_product_quadrature)
      {
        pt.put ("settings.problem.quadrature.azimuthal_order", m_order[0]);
        pt.put ("settings.problem.quadrature.polar_order", m_order[1]);
      }
      else
        pt.put ("settings.problem.quadrature.sn", m_order[0]);
    }
  }

  // InputAlgebra --------------------------------------------------------------

  InputAlgebra::InputAlgebra ()
      : m_matrix_free (true),
        m_use_fission_density (true),
        m_eig_solver ("PI"),
        m_eig_tol (1.e-6),
        m_eig_max_it (1000),
        m_mg_solver ("GS"),
        m_mg_tol (1.e-6),
        m_mg_max_it (1000),
        m_inner_solver ("Krylov"),
        m_inner_tol (m_eig_tol / 10.),
        m_inner_max_it (1000)
  {
    m_n_cv = 3;
    m_n_eigenvalues = 1;
    m_use_power_iteration = true;
  }

  void
  InputAlgebra::set (bool matrix_free,
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
    m_matrix_free = matrix_free;

    m_use_fission_density = use_fission_density;

    m_eig_solver = eig_solver;
    m_eig_tol = eig_tol;
    m_eig_max_it = eig_max_it;

    m_mg_solver = mg_solver;
    m_mg_tol = mg_tol;
    m_mg_max_it = mg_max_it;

    m_inner_solver = inner_solver;
    m_inner_tol = inner_tol;
    m_inner_max_it = inner_max_it;
  }

  void
  InputAlgebra::read_xml (const boost::property_tree::ptree & pt,
                      const std::string & /* path */,
                      const std::string & /* filename_no_ext */)
  {
    m_matrix_free = pt.get<bool> ("settings.algebra.matrix_free");

    m_eig_solver = pt.get<std::string> (
        "settings.algebra.eig_solver.<xmlattr>.type", "PI");

    if (m_eig_solver.compare (std::string ("PI")) == 0)
    {
      m_use_power_iteration = true;
    }
    else if (m_eig_solver.length () > 6
        and
        (m_eig_solver.compare (0, 6, "krylov") == 0 or
         m_eig_solver.compare (0, 6, "Krylov") == 0))
    {
      m_use_power_iteration = false;

      std::size_t pos1 = m_eig_solver.find_first_of(":");
      std::size_t pos2 = m_eig_solver.find_last_of(":");

      if (m_eig_solver.substr (7).npos == 0)
      {
        m_n_eigenvalues = 1;
        m_n_cv = 2 * m_n_eigenvalues + 1;
      }
      else
      {
        m_n_eigenvalues = string_to_num<unsigned int> (m_eig_solver.substr(pos1+1, pos2-pos1-1));
        m_n_cv = string_to_num<unsigned int> (m_eig_solver.substr(pos2+1));
      }
    }
    else
    {
      // TODO Throw an error
      std::exit (1);
    }

    std::string bin = pt.get<std::string> (
        "settings.algebra.eig_solver.<xmlattr>.form", "FD");
    if (bin.compare (std::string ("FD")) == 0)
      m_use_fission_density = true;
    else
      m_use_fission_density = false;

    m_eig_tol = pt.get<double> ("settings.algebra.eig_solver.tol", 1.e-6);
    m_eig_max_it = pt.get<unsigned int> ("settings.algebra.eig_solver.max_it",
        1000);

    m_inner_solver = pt.get<std::string> (
        "settings.algebra.inner_solver.<xmlattr>.type", "Krylov");
    m_inner_tol = pt.get<double> ("settings.algebra.inner_solver.tol",
        m_eig_tol / 10.);
    m_inner_max_it = pt.get<unsigned int> ("settings.algebra.inner_solver.max_it",
        1000);

    m_mg_solver = pt.get<std::string> (
        "settings.algebra.mg_solver.<xmlattr>.type", "GS");
    m_mg_tol = pt.get<double> ("settings.algebra.mg_solver.tol",
        m_eig_tol / 10.);
    m_mg_max_it = pt.get<unsigned int> ("settings.algebra.mg_solver.max_it",
        100);

  }

  void
  InputAlgebra::write (boost::property_tree::ptree & pt) const
  {

    pt.add ("settings.<xmlcomment>", "algebra");
    pt.add ("settings.algebra", "");

    pt.put ("settings.algebra.matrix_free", m_matrix_free);

    pt.put ("settings.algebra.eig_solver.<xmlattr>.type", m_eig_solver);
    pt.put ("settings.algebra.eig_solver.tol", m_eig_tol);
    pt.put ("settings.algebra.eig_solver.max_it", m_eig_max_it);

    pt.put ("settings.algebra.inner_solver.<xmlattr>.type", m_inner_solver);
    pt.put ("settings.algebra.inner_solver.tol", m_inner_tol);
    pt.put ("settings.algebra.inner_solver.max_it", m_inner_max_it);
  }

  // InputFESettings -----------------------------------------------------------

  InputFESettings::InputFESettings ()
      : m_degree (0),
        m_n_ref (0)
  {
  }

  void
  InputFESettings::set (unsigned int degree, unsigned int n_ref)
  {
    m_degree = degree;
    m_n_ref = n_ref;
  }

  void
  InputFESettings::read_xml (const boost::property_tree::ptree & pt,
                         const std::string & /* path */,
                         const std::string & /* filename_no_ext */)
  {
    m_degree = pt.get<unsigned int> ("settings.fe_settings.<xmlattr>.degree", 0);
    m_n_ref = pt.get<unsigned int> ("settings.fe_settings.<xmlattr>.n_ref", 0);
  }

  void
  InputFESettings::write (boost::property_tree::ptree & pt) const
  {
    pt.add ("settings.<xmlcomment>", "fe_settings");
    pt.put ("settings.fe_settings.<xmlattr>.degree", m_degree);
    pt.put ("settings.fe_settings.<xmlattr>.n_ref", m_n_ref);
  }

  // InputOutput ---------------------------------------------------------------

  InputOutput::InputOutput ()
      : homogenize (false),
        print_mesh (false),
        vtk_power (false),
        vtk_flux (false),
        vtk_mat (false),
        ass_bcs (false)
  {
  }

  void
  InputOutput::set (const std::vector<bool> & options)
  {
    homogenize = options[0];
    print_mesh = options[1];
    vtk_power = options[2];
    vtk_flux = options[3];
    vtk_mat = options[4];
    ass_bcs = options[5];
  }

  void
  InputOutput::read_xml (const boost::property_tree::ptree & pt,
                     const std::string & /* path */,
                     const std::string & /* filename_no_ext */)
  {
    homogenize = pt.get<bool> ("settings.output.homogenize", false);
    print_mesh = pt.get<bool> ("settings.output.print_mesh", false);
    vtk_power = pt.get<bool> ("settings.output.vtk_power", false);
    vtk_flux = pt.get<bool> ("settings.output.vtk_flux", false);
    vtk_mat = pt.get<bool> ("settings.output.vtk_mat", false);
    ass_bcs = pt.get<bool> ("settings.output.ass_bcs", false);
  }

  void
  InputOutput::write (boost::property_tree::ptree & pt) const
  {
    pt.add ("settings.<xmlcomment>", "output");
    pt.add ("settings.output", "");

    pt.put ("settings.output.homogenize", homogenize);
    pt.put ("settings.output.print_mesh", print_mesh);
    pt.put ("settings.output.vtk_power", vtk_power);
    pt.put ("settings.output.vtk_flux", vtk_flux);
    pt.put ("settings.output.vtk_mat", vtk_mat);
    pt.put ("settings.output.ass_bcs", ass_bcs);
  }

  // InputSettings -------------------------------------------------------------

  InputSettings::InputSettings ()
        : settings_ext(".settings.xml"),
          files (),
        problem (),
        algebra (),
        fe_settings (),
        output ()
  {
  }

  InputSettings::InputSettings (const std::string & path_plus_file_plus_ext_)
      : settings_ext(std::string(".settings.xml")),
        path_plus_file_plus_ext (path_plus_file_plus_ext_),
        path (Forest::get_path (path_plus_file_plus_ext)),
        problem_name (Forest::get_file_plus_ext (path_plus_file_plus_ext)),
        files (),
        problem (),
        algebra (),
        fe_settings (),
        output ()
  {
    this->load_xml (path, problem_name);
    // this->save (path, problem_name);
  }

  void
  InputSettings::load_xml (const std::string &path_,
                       const std::string &problem_name_)
  {
    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    // Load XML file and put its contents in property tree. 
    // No namespace qualification is needed, because of Koenig 
    // lookup on the second argument. If reading fails, exception
    // is thrown.
    read_xml (path_ + problem_name_ + settings_ext, pt,
        boost::property_tree::xml_parser::trim_whitespace);


    files.read_xml (pt, path_, problem_name_);
    //coarsening_bcs.read_xml (pt, path_, problem_name_);
    output.read_xml (pt, path_, problem_name_);
    problem.read_xml (pt, path_, problem_name_);
    algebra.read_xml (pt, path_, problem_name_);
    fe_settings.read_xml (pt, path_, problem_name_);
  }

  void
  InputSettings::save_xml (const std::string &path_,
                       const std::string &problem_name_)
  {
    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    pt.add ("settings", "");
    pt.add ("settings.<xmlcomment>",
        " \n \
      Here we define the different materials associated to the \n \
      material id specified when defining the geometry");

    // the rest
    files.write (pt);
    files.write (pt);
    output.write (pt);
    problem.write (pt);
    algebra.write (pt);
    fe_settings.write (pt);

    // Write property tree to XML file (\t is a tabulator)
    write_xml (path_ + "out_" + problem_name_ + settings_ext, pt,
        std::locale (), xml_writer_settings('\t', 1));
  }

} // end of namespace Forest
