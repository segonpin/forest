/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   output/forest_post_processing.cc
 * @brief  Implementation of class template PostProcessing
 */

#include "output/forest_post_processing.h"

#include "utils/forest_utils_base.h"
#include "input/input.h"
#include "geometry/prob_geom.h"
#include "neutronics/state.h"
#include "algebra/equation.h"

#include "utils/forest_utils_logstream.h" // for log_print_text
#include "utils/forest_utils_memory.h"   // for StatsMemory
#include "utils/forest_utils_timer.h"    // for StatsTimer

#include "deal.II/base/timer.h"          // for Timer
#include "deal.II/dofs/dof_handler.h"    // for DoFHandler
#include "deal.II/base/logstream.h"
#include "deal.II/base/quadrature_lib.h" // for QGauss
#include "deal.II/grid/tria.h"
#include "deal.II/numerics/data_out.h"
#include "deal.II/base/exceptions.h"     // for Assert
#include "deal.II/lac/vector.h"          // for Vector
#include "deal.II/fe/fe_values.h"        // for FEValues
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <fstream>                      // for basic_filebuf<>::int_type, etc
#include <sstream>                      // for basic_stringbuf<>::int_type, etc
#include <cmath>
#include <iostream>
#include <locale>
#include <string>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  PostProcessing<dim>::PostProcessing (const State<dim> & state)
      : mp_state (state),
        data (mp_state.mp_data),
        geom (mp_state.mp_geom),
        fe (mp_state.get_fe ()),
        dof_handler (mp_state.get_dof_handler ()),
        m_power_dofs (mp_state.get_n_elem ()),
        m_scaling (1.0),
        m_verbose (false)
  {
    // start the timer for initializing the geometry object.
    Timer timer;
    timer.restart (); // We restart the timer

    // message
    if (m_verbose)
      log_print_text ("PostProcessing", 1);

    build_flux_power_data ();
    print_xml_data (data.mp_settings.get_file_out ());
    print_image_data ();

    // Stop the timer and report the time
    timer.stop ();
    StatsTimer::instance ().add ("PostProcessing", timer.wall_time ());
    // Reprot the memory
    //StatsMemory::instance ().add ("PostProcessing", memory_consumption ());
  }

  /**
   * @todo add weights to consider the different importance of every
   * cell when adding for the pin, and the importance of every pin when adding
   * for the assembly.
   * @todo find a normalization for the power, other for the fluxes,
   * and write which one is used and how can it be used to move from one to the
   * other (keep the scaling factors?)
   */
  template <int dim>
  void
  PostProcessing<dim>::build_flux_power_data ()
  {
    if (m_verbose)
      std::cout << "start build_pin_power" << std::endl;

    // We need this constants
    const unsigned int n_groups = data.mp_mat.get_n_groups ();
    const unsigned int n_dofs = dof_handler.n_dofs();
    const unsigned int n_cells = geom.m_tria.n_active_cells ();
    const unsigned int n_pins = geom.get_n_pins ();
    const unsigned int n_ass = geom.get_n_assemblies ();

    // Initialize all that  is needed to iterate over dofs and cells
    const unsigned int fe_degree = geom.mp_data.mp_settings.get_fe_degree ();
    QGauss<dim> quadrature_formula (fe_degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula,
        update_values | update_quadrature_points | update_volume_elements
        | update_JxW_values);
    const unsigned int n_q_points = quadrature_formula.size ();
    std::vector<unsigned int> global_dof_ind (fe.dofs_per_cell);

    // Now, let's multiply the scalar flux per group per the
    // corresponding fission cross-section, and add them to get
    // the power distribution.
    for (typename DoFHandler<dim, dim>::active_cell_iterator cell_it =
        dof_handler.begin_active (); cell_it != dof_handler.end (); ++cell_it)
    {
      fe_values.reinit (cell_it);
      cell_it->get_dof_indices (global_dof_ind);
      const unsigned int mat = cell_it->material_id ();
      std::vector<double> nusigfission (n_groups, 0.);
      if (data.mp_mat.m_xs.at (mat).exist_sigf)
      {
        nusigfission = data.mp_mat.m_xs.at (mat).m_sigf;
      }
      else if (data.mp_mat.m_xs.at (mat).exist_nusigf)
      {
        nusigfission = data.mp_mat.m_xs.at (mat).m_nusigf;
      }
      else
      {
        nusigfission.resize (n_groups, 0.);
      }
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
        {
          m_power_dofs[global_dof_ind[i]] += nusigfission[g]
              * mp_state.m_scalar_flux.m_data[g].block (0)[global_dof_ind[i]];
        }
      }
    }

    // We go from dof to cell -------------------------------------------------

    // how much is the weight of each pin?
    m_weight_cell.resize (n_cells, 0.0);
    for (typename DoFHandler<dim, dim>::active_cell_iterator cell_it =
        dof_handler.begin_active (); cell_it != dof_handler.end (); ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      m_weight_cell[cell] = cell_it->measure ();
    }
    // Now we get the integral of the flux and power over each cell
    m_power_cell.resize (n_cells, 0.0);
    m_flux_cell.resize (n_groups, std::vector<double> (n_cells, 0.0));
    for (typename DoFHandler<dim, dim>::active_cell_iterator cell_it =
        dof_handler.begin_active (); cell_it != dof_handler.end (); ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      fe_values.reinit (cell_it);
      cell_it->get_dof_indices (global_dof_ind);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        // Calculate the weight of this degree of freedom
        double weight_dof_i = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          weight_dof_i += fe_values.shape_value (i, q) * fe_values.JxW (q);
        }
        // add the corresponding contribution
        m_power_cell[cell] += m_power_dofs[global_dof_ind[i]] * weight_dof_i;
        for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
        {
          m_flux_cell[g][cell] += weight_dof_i
              * mp_state.m_scalar_flux.m_data[g].block (0)[global_dof_ind[i]];
        }
      }
    }

    // We go from cell to pin -------------------------------------------------

    // how much is the weight of each pin?
    m_weight_pin.resize (n_pins, 0);
    for (typename DoFHandler<dim, dim>::active_cell_iterator cell_it =
        dof_handler.begin_active (); cell_it != dof_handler.end (); ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = geom.get_cell2pin (cell);
      m_weight_pin[pin] += m_weight_cell[cell];
    }
    // Now we get the flux and power over each pin
    m_power_pin.resize (n_pins, 0);
    m_flux_pin.resize (n_groups, std::vector<double> (n_pins, 0.0));
    for (typename DoFHandler<dim, dim>::active_cell_iterator cell_it =
        dof_handler.begin_active (); cell_it != dof_handler.end (); ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = geom.get_cell2pin (cell);
      m_power_pin[pin] += m_power_cell[cell];
      for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
      {
        m_flux_pin[g][pin] += m_flux_cell[g][cell];
      }
    }

    // We go from pin to assembly ---------------------------------------------
    ///@bug if m_flux_pin is already multiplied by the length of the interval
    /// "weight", then it should not be multplied again by this when averaging
    /// for m_flux_ass . I should check that this is ok.

    // how much is the weight of each assembly?
    m_weight_ass.resize (n_ass, 0);
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    {
      const unsigned int ass = geom.get_pin2assembly (pin);
      m_weight_ass[ass] += m_weight_pin[pin];
    }
    // Now we get the flux and power over each assembly
    m_power_ass.resize (n_ass, 0);
    m_flux_ass.resize (n_groups, std::vector<double> (n_ass, 0.0));
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    {
      const unsigned int ass = geom.get_pin2assembly (pin);
      m_power_ass[ass] += m_power_pin[pin];
      for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
      {
        m_flux_ass[g][ass] += m_flux_pin[g][pin];
      }
    }

    // We go from assembly to total---------------------------------------------
    double power_total = 0.0;
    double weight_total = 0.0;
    for (unsigned int i = 0; i < m_power_ass.size (); ++i)
    {
      weight_total += m_weight_ass[i];
      power_total += m_power_ass[i];
    }

    // ------------------------------------------------------------------------
    // normalization ----------------------------------------------------------
    // pow*w' / sum(w)
    m_scaling = 1. / (power_total/weight_total);


    // Normalization (average power per dofs is equal to one) ------------------
    m_norm_power_dofs.reinit (n_dofs, 0.0);
    m_norm_flux_dofs.resize (n_groups, m_norm_power_dofs);
    for (unsigned int i = 0; i < m_power_dofs.size (); ++i)
    {
      m_norm_power_dofs[i] = m_power_dofs[i]*m_scaling;
    }
    // Normalization (average flux per ass is equal to one (for fast group))
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        m_norm_flux_dofs[g][i] = mp_state.m_scalar_flux.m_data[g].block (0)[i]
            * m_scaling;
      }
    }

    // Normalization (average power per cell is equal to one) ------------------
    m_norm_power_cell.resize (n_cells, 0.0);
    m_norm_flux_cell.resize (n_groups, std::vector<double> (n_cells, 0.0));
    for (unsigned int i = 0; i < m_power_cell.size (); ++i)
    {
      m_norm_power_cell[i] = m_power_cell[i]/m_weight_cell[i]*m_scaling;
    }
    // Normalization (average flux per ass is equal to one (for fast group))
    for (unsigned int g = 0; g < m_flux_cell.size (); ++g)
    {
      for (unsigned int i = 0; i < m_flux_cell[g].size (); ++i)
      {
        m_norm_flux_cell[g][i] = m_flux_cell[g][i] / m_weight_cell[i]
            * m_scaling;
      }
    }

    // Normalization (average power per pin is equal to one) ------------------
    m_norm_power_pin.resize (n_pins, 0);
    m_norm_flux_pin.resize (n_groups, std::vector<double> (n_pins, 0.0));
    for (unsigned int i = 0; i < m_power_pin.size (); ++i)
    {
      m_norm_power_pin[i] = m_power_pin[i]/m_weight_pin[i]*m_scaling;
    }
    // Normalization (average flux per ass is equal to one (for fast group))
    for (unsigned int g = 0; g < m_flux_pin.size (); ++g)
    {
      for (unsigned int i = 0; i < m_flux_pin[g].size (); ++i)
      {
        m_norm_flux_pin[g][i] = m_flux_pin[g][i]/m_weight_pin[i]*m_scaling;
      }
    }

    // Normalization (average power per ass is equal to one) ------------------
    m_norm_power_ass.resize (n_pins, 0);
    m_norm_flux_ass.resize (n_groups, std::vector<double> (n_ass, 0.0));
    for (unsigned int i = 0; i < m_power_ass.size (); ++i)
    {
      m_norm_power_ass[i] = m_power_ass[i]/m_weight_ass[i]*m_scaling;
    }
    // Normalization (average flux per ass is equal to one (for fast group))
    for (unsigned int g = 0; g < m_flux_ass.size (); ++g)
    {
      for (unsigned int i = 0; i < m_flux_ass[g].size (); ++i)
      {
        m_norm_flux_ass[g][i] = m_flux_ass[g][i]/m_weight_ass[i]*m_scaling;
      }
    }

  }

  template <int dim>
  double
  PostProcessing<dim>::calculate_scaling (const std::vector<double> & vec)
  {
    // Normalization (average power per pin is equal to one)
    double vec_sum = 0;
    double vec_counter = 0;
    for (unsigned int i = 0; i < vec.size (); ++i)
    {
      if (std::abs (vec[i]) > 1.e-10)
      {
        vec_sum += vec[i];
        vec_counter += 1.0;
      }
    }
    double scal = vec_counter / vec_sum;
    return scal;
  }

  template <int dim>
  void
  PostProcessing<dim>::print_pin_data (std::string & bin,
                                       const std::vector<double> & quantity,
                                       const std::string & indent,
                                       const unsigned int n_indent)
  {
    if (m_verbose)
      std::cout << "start print_pin_data" << std::endl;

    // Changing the container
    // m_power_pin.resize (geom.tria.n_active_cells (), 0);
    std::vector<std::vector<double> > bin_power_pin;
    std::vector<double> tmp_power_pin;
    for (unsigned int j = 0; j < geom.map_pin_index.size (); ++j)
    {
      for (unsigned int i = 0; i < geom.map_pin_index[j].size (); ++i)
      {
        tmp_power_pin.push_back (quantity[geom.map_pin_index[j][i]]);
      }
      bin_power_pin.push_back (tmp_power_pin);
      tmp_power_pin.clear ();
    }

    bin.clear ();
    Forest::vector_to_string (bin_power_pin, bin, indent, n_indent);
  }

  template <int dim>
  void
  PostProcessing<dim>::print_ass_data (std::string & bin,
                                       const std::vector<double> & quantity,
                                       const std::string &indent,
                                       const unsigned int n_indent)
  {
    if (m_verbose)
      std::cout << "start print_ass_data" << std::endl;

    // Changing the container
    // m_power_pin.resize (geom.tria.n_active_cells (), 0);
    std::vector<std::vector<double> > bin_power_ass;
    std::vector<double> tmp_power_ass;
    for (unsigned int j = 0; j < geom.m_core.m_indices[0].size (); ++j)
    {
      for (unsigned int i = 0; i < geom.m_core.m_indices[0][j].size (); ++i)
      {
        tmp_power_ass.push_back (quantity[geom.m_core.m_indices[0][j][i]]);
      }
      bin_power_ass.push_back (tmp_power_ass);
      tmp_power_ass.clear ();
    }

    bin.clear ();
    Forest::vector_to_string (bin_power_ass, bin, indent, n_indent);
  }

  template <int dim>
  void
  PostProcessing<dim>::print_xml_data (const std::string &filename)
  {
    if (m_verbose)
      std::cout << "start print_xml_data" << std::endl;

    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    std::string bin; // bin is a container to save temporal data
    std::string indent ("\t");
    const unsigned int n_indent = 3;
    std::string field_root;
    std::string field;

    // We print the keff ------------------------------------------------------
    pt.add ("results", "");
    pt.add ("results.<xmlcomment>", "K-effective for the reactor");
    pt.put ("results.keff", mp_state.get_keff ());

    // We print information about the pins ------------------------------------
    pt.add ("results.<xmlcomment>", "Retults at pin-level");
    field_root = "results.pin";
    pt.add (field_root, "");
    // pt.put (field_root + ".<xmlattr>.sym", false);

    field = field_root + ".weight";
    this->print_pin_data (bin, m_weight_pin, indent, n_indent);
    pt.put (field, bin);
    bin.clear ();

    field = field_root + ".power";
    this->print_pin_data (bin, m_norm_power_pin, indent, n_indent);
    pt.put (field, bin);
    bin.clear ();

    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      field = field_root + ".flux_" + Forest::num_to_string (g);
      this->print_pin_data (bin, m_norm_flux_pin[g], indent, n_indent);
      pt.put (field, bin);
      bin.clear ();
    }

    // We print information about the assemblies ------------------------------
    pt.add ("results.<xmlcomment>", "Retults at assembly-level");
    field_root = "results.ass";
    pt.add (field_root, "");
    // pt.put (field_root + ".<xmlattr>.sym", false);

    field = field_root + ".weight";
    this->print_ass_data (bin, m_weight_ass, indent, n_indent);
    pt.put (field, bin);
    bin.clear ();

    field = field_root + ".power";
    this->print_ass_data (bin, m_norm_power_ass, indent, n_indent);
    pt.put (field, bin);
    bin.clear ();

    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      field = field_root + ".flux_" + Forest::num_to_string (g);
      this->print_ass_data (bin, m_norm_flux_ass[g], indent, n_indent);
      pt.put (field, bin);
      bin.clear ();
    }

    // Write property tree to XML file (\t is a tabulator) --------------------
#if BOOST_VERSION < 105500
    boost::property_tree::xml_writer_settings<char> settings ('\t', 1);
#else
    boost::property_tree::xml_writer_settings<std::string> settings ('\t', 1);
#endif
    write_xml (filename, pt, std::locale (), settings);
  }

  template <int dim>
  void
  PostProcessing<dim>::print_vtk () const
  {
    if (m_verbose)
      std::cout << "start print_vtk" << std::endl;

    // Preparing the vtk file
    std::string filename = data.mp_settings.get_file_vtk ();
    Forest::log_print_text ("Writing solution to <" + filename + ">", 1);
    std::ofstream vtk_output (filename.c_str ());
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);

    // We add the fission rate
    data_out.add_data_vector (m_norm_power_dofs, "power",
        DataOut<dim>::type_dof_data);

    // Here we add the scalar flux
    for (unsigned int g = 0; g < mp_state.get_n_groups (); ++g)
    {
      std::string field_name = "Phi_" + Forest::num_to_string (g);
      data_out.add_data_vector (m_norm_flux_dofs[g],
          field_name.c_str (), DataOut<dim>::type_dof_data);
    }

    typename Triangulation<dim, dim>::active_cell_iterator cell;

    // Here we add the materials for every cell
    Vector<double> mat_id_cells (geom.m_tria.n_active_cells ());
    cell = geom.m_tria.begin_active ();
    for (unsigned int ind = 0; cell != geom.m_tria.end (); ++cell, ++ind)
      mat_id_cells[ind] = static_cast<double> (cell->material_id ());
    data_out.add_data_vector (mat_id_cells, "material",
        DataOut<dim>::type_cell_data);

    // Here we add the user index for every cell
    Vector<double> user_index (geom.m_tria.n_active_cells ());
    cell = geom.m_tria.begin_active ();
    for (unsigned int ind = 0; cell != geom.m_tria.end (); ++cell, ++ind)
      user_index[ind] = static_cast<double> (cell->user_index ());
    data_out.add_data_vector (user_index, "user_index",
        DataOut<dim>::type_cell_data);

    // Here we add the pin index for every cell
    Vector<double> pin_index (geom.m_tria.n_active_cells ());
    Vector<double> ass_index (geom.m_tria.n_active_cells ());
    cell = geom.m_tria.begin_active ();
    for (unsigned int ind = 0; cell != geom.m_tria.end (); ++cell, ++ind)
    {
      const unsigned int cell_ind = cell->user_index ();
      const unsigned int pin_ind = mp_state.mp_geom.get_cell2pin (cell_ind);
      pin_index[ind] = static_cast<double> (pin_ind);
      const unsigned int ass_ind = mp_state.mp_geom.get_pin2assembly (pin_ind);
      ass_index[ind] = static_cast<double> (ass_ind);
    }
    data_out.add_data_vector (pin_index, "pin_index",
        DataOut<dim>::type_cell_data);
    data_out.add_data_vector (ass_index, "assembly_index",
        DataOut<dim>::type_cell_data);

    // We refine the mesh depending on the fe_degree, so that we do
    const unsigned int fe_degree = geom.mp_data.mp_settings.get_fe_degree ();
    const unsigned int n_patches =
        fe_degree == 0 ? 0 : std::pow (2, fe_degree - 1);
    data_out.build_patches (n_patches);

    // We print the data to the file
    data_out.write_vtk (vtk_output);
  }

  template <int dim>
  void
  PostProcessing<dim>::print_eps () const
  {
    if (m_verbose)
      std::cout << "start print_eps" << std::endl;

    /*
     //Write the grid in eps format if we are in more than 1d
     if (dim > 1 && false)
     {
     std::string filename = "grid";
     filename += ".eps";
     deallog << "Writing grid to <" << filename << ">" << std::endl;
     std::ofstream eps_output (filename.c_str ());
     GridOut grid_out;
     grid_out.write_eps (geom.tria, eps_output);
     }
     */

  }

  template <int dim>
  void
  PostProcessing<dim>::print_gpl () const
  {
    if (m_verbose)
      std::cout << "start print_gpl" << std::endl;

    /*
     //Output of the solution in gnuplot format.
     std::string filename = "sol";
     filename += ".gnuplot";
     deallog << "Writing solution to <" << filename << ">" << std::endl;
     std::ofstream gnuplot_output (filename.c_str ());
     DataOut<dim> data_out;
     data_out.attach_dof_handler (dof_handler);
     data_out.add_data_vector (solution.block (0), "u");
     data_out.build_patches ();
     data_out.write_gnuplot (gnuplot_output);
     */

  }

  template <int dim>
  void
  PostProcessing<dim>::print_image_data () const
  {
    if (m_verbose)
      std::cout << "start output_results" << std::endl;

    // Printing to the vtk file
    this->print_vtk ();
    // Printing to the eps file
    this->print_eps ();
    // Printing to the gnuplot file
    this->print_gpl ();
  }

  template class PostProcessing<1> ;
  template class PostProcessing<2> ;
  template class PostProcessing<3> ;

} // end of namespace Forest
