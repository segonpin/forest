/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/homogenization.cc
 * @brief  Implementation of class template Homogenization
 */

#include "neutronics/homogenization.h"

#include "utils/forest_utils_base.h"
#include "input/input.h"
#include "input/input_settings.h"
#include "input/input_mat.h"
#include "geometry/prob_geom.h"
#include "neutronics/state.h"
#include "neutronics/extract_bv.h"      // for ExtractBv

#include "deal.II/base/quadrature_lib.h"  // for QGauss
#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_values.h"       // for FEValues
#include "deal.II/fe/fe_update_flags.h"  // for operator|, etc
#include "deal.II/base/exceptions.h"    // for Assert and ExcMessage

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <string>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  Homogenization<dim>::Homogenization (const State<dim> & state/*,
   const ExtractBV<dim> & extract_bv*/)
      : mp_state (state),
        /*mp_extract_bv (extract_bv),*/
        data (mp_state.mp_data),
        geom (mp_state.mp_geom),
        fe (mp_state.get_fe ()),
        dof_handler (mp_state.get_dof_handler ())
  {
    /* First we calculate the flux weight per cell  */
    cell_weights ();

    /* First we homogenize */
    pin_homogenization ();

    /* First we homogenize */
    assembly_homogenization ();

    /*
     // We make this members of extract_bv visible in this function.
     const stdvd3 & phi0_pin = mp_extract_bv.m_phi0_pin;
     const stdvd3 & phi2_pin = mp_extract_bv.m_phi2_pin;
     const stdvd3 & current_pin = mp_extract_bv.m_j_pin;
     // now we print the values
     std::string f_pin_xs (
     data.mp_settings.get_path_plus_file () + ".pin.mat.xml");
     print_homogenized_xs (f_pin_xs, pin_xs, phi0_pin, phi2_pin, current_pin);

     // We make this members of extract_bv visible in this function.
     const stdvd3 & phi0_ass = mp_extract_bv.m_phi0_ass;
     const stdvd3 & phi2_ass = mp_extract_bv.m_phi2_ass;
     const stdvd3 & current_ass = mp_extract_bv.m_j_ass;
     // now we print the values
     std::string f_ass_xs (
     data.mp_settings.get_path_plus_file () + ".ass.mat.xml");
     print_homogenized_xs (f_ass_xs, ass_xs, phi0_ass, phi2_ass, current_ass);*/

    /**
     @todo Use ExtractBV to do this
     @code{.cpp}
     if (data.mp_settings.get_print_ass_bcs_flag ())
     {
     state.print_bv (data.mp_settings.get_path_plus_file_plus_ext ());
     }
     @endcode
     */
  }

  template <int dim>
  void
  Homogenization<dim>::cell_weights ()
  {
    /* We need some constants */
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_cells = geom.m_tria.n_active_cells ();

    /* Initialize all that  is needed to iterate over dofs and cells */
    unsigned int fe_degree = data.mp_settings.get_fe_degree ();
    QGauss<dim> quad_formula (fe_degree + 1);
    FEValues<dim> fe_values (fe, quad_formula,
        update_values | update_quadrature_points | update_volume_elements
        | update_JxW_values);

    const unsigned int n_q = quad_formula.size ();
    std::vector<unsigned int> global_dof_ind (fe.dofs_per_cell);

    /* Assuming the cross sections is constant per cell, we need the
     * integral of the flux per cell, that will be used for the weighting
     * of the cross sections */
    cell_w.resize (n_cells, std::vector<double> (n_groups, 0.));

    /* For shorter access to the scalar flux */
    const SnVector & flux = mp_state.m_scalar_flux;

    /* Run over all the cells to calculate each weight */
    c_iter cell_it = dof_handler.begin_active (), endc = dof_handler.end ();
    for (; cell_it != endc; ++cell_it)
    {
      fe_values.reinit (cell_it);
      cell_it->get_dof_indices (global_dof_ind);
      const unsigned int cell = cell_it->user_index ();

      /* for every degree of freedom in this cell... */
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        /* ...we calculate the integral of the associated base function, */
        double integral = 0;
        for (unsigned int q = 0; q < n_q; ++q)
          integral += fe_values.shape_value (i, q) * fe_values.JxW (q);

        /* and then we use it for the density of every group. */
        for (unsigned int g = 0; g < n_groups; ++g)
          cell_w[cell][g] += flux.m_data[g][global_dof_ind[i]] * integral;
      }
    }
  }

  template <int dim>
  void
  Homogenization<dim>::pin_homogenization ()
  {
    /* We need some constants */
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_pins = geom.get_n_pins ();

    /* Container for the number of pins times the number of groups */
    pin_w.resize (n_pins, std::vector<double> (n_groups, 0.));
    // Container for the homogenized cross section
    pin_xs.resize (n_pins);

    // Calculate the pin weight of the flux
    c_iter cell_it, endc = dof_handler.end ();
    for (cell_it = dof_handler.begin_active (); cell_it != endc; ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = geom.get_cell2pin (cell);
      for (unsigned int g = 0; g < n_groups; ++g)
        pin_w[pin][g] += cell_w[cell][g];
    }

    // Initialize to zero the container for the homogenized cross sections
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    {
      this->init_xs (pin_xs[pin], n_groups);
    }

    // Run over all the cells to calculate each weight
    for (cell_it = dof_handler.begin_active (); cell_it != endc; ++cell_it)
    {
      const unsigned int mat = cell_it->material_id ();
      const unsigned int cell = cell_it->user_index ();
      const unsigned int pin = geom.get_cell2pin (cell);
      this->add_cell_contribution (pin_xs[pin], data.mp_mat.m_xs.at (mat), cell_w[cell]);
    }

    // Normalize dividing by the integral of the scalar flux over the pin
    for (unsigned int pin = 0; pin < n_pins; ++pin)
    {
      this->normalize (pin_xs[pin], pin_w[pin]);
    }
  }

  template <int dim>
  void
  Homogenization<dim>::assembly_homogenization ()
  {
    // We need some constants
    const unsigned int n_groups = mp_state.get_n_groups ();
    const unsigned int n_ass = geom.get_n_assemblies ();

    // Container for the number of pins times the number of groups
    ass_w.resize (n_ass, std::vector<double> (n_groups, 0.));
    // Container for the homogenized cross section
    ass_xs.resize (n_ass);

    // Calculate the assembly weight of the flux
    c_iter cell_it = dof_handler.begin_active (), endc = dof_handler.end ();
    for (cell_it = dof_handler.begin_active (); cell_it != endc; ++cell_it)
    {
      const unsigned int cell = cell_it->user_index ();
      const unsigned int ass = geom.get_cell2assembly (cell);
      for (unsigned int g = 0; g < n_groups; ++g)
        ass_w[ass][g] += cell_w[cell][g];
    }

    // Initialize to zero the container for the homogenized cross sections
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    {
      this->init_xs (ass_xs[ass], n_groups);
    }

    // Run over all the cells to calculate each weight
    for (cell_it = dof_handler.begin_active (); cell_it != endc; ++cell_it)
    {
      const unsigned int mat = cell_it->material_id ();
      const unsigned int cell = cell_it->user_index ();
      const unsigned int ass = geom.get_cell2assembly (cell);
      this->add_cell_contribution (ass_xs[ass], data.mp_mat.m_xs.at (mat), cell_w[cell]);
    }

    // Normalize dividing by the integral of the scalar flux over the pin
    for (unsigned int ass = 0; ass < n_ass; ++ass)
    {
      this->normalize (ass_xs[ass], ass_w[ass]);
    }
  }

  template <int dim>
  void
  Homogenization<dim>::init_xs (XS_single & dst,
                                const unsigned int n_groups)
  {
    dst.n_groups = n_groups;
    dst.m_sigmat.resize (n_groups, 0.);
    dst.m_sigmas.resize (n_groups, std::vector<double> (n_groups, 0.));
    dst.m_nusigf.resize (n_groups, 0.);
    dst.m_sigf.resize (n_groups, 0.);
    dst.m_chi.resize (n_groups, 0.);
  }

  template <int dim>
  void
  Homogenization<dim>::add_cell_contribution (XS_single & dst,
                                              const XS_single & src,
                                              const std::vector<double> & weight)
  {
    // Verify that the data exists
    Assert(src.exist_sigmat, ExcMessage("sigmat not initialized"));
    Assert(src.exist_sigmas, ExcMessage("sigmas not initialized"));
    Assert(src.exist_chi, ExcMessage("chi not initialized"));
    Assert(src.exist_nusigf, ExcMessage("nusigf not initialized"));

    // Add the cell contribution to the cross section for this pin
    const unsigned int n_groups = src.n_groups;
    for (unsigned int g = 0; g < n_groups; ++g)
      dst.m_sigmat[g] += weight[g] * src.m_sigmat[g];
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int h = 0; h < n_groups; ++h)
        dst.m_sigmas[g][h] += weight[h] * src.m_sigmas[g][h];
    for (unsigned int g = 0; g < n_groups; ++g)
      dst.m_nusigf[g] += weight[g] * src.m_nusigf[g];
    if (src.exist_sigf)
      for (unsigned int g = 0; g < n_groups; ++g)
        dst.m_sigf[g] += weight[g]* src.m_sigf[g];
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int h = 0; h < n_groups; ++h)
        dst.m_chi[g] += src.m_chi[g] * src.m_nusigf[h] * weight[h];
  }

  template <int dim>
  void
  Homogenization<dim>::normalize (XS_single & dst,
                                  const std::vector<double> & weight)
  {
    const unsigned int n_groups = dst.n_groups;

    dst.exist_sigmat = true;
    dst.exist_sigmas = true;
    dst.exist_chi = true;
    dst.exist_nusigf = true;
    dst.exist_sigf = true;

    for (unsigned int g = 0; g < n_groups; ++g)
      dst.m_sigmat[g] = dst.m_sigmat[g] / weight[g];
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int h = 0; h < n_groups; ++h)
        dst.m_sigmas[g][h] = dst.m_sigmas[g][h] / weight[h];
    for (unsigned int g = 0; g < n_groups; ++g)
      dst.m_nusigf[g] = dst.m_nusigf[g] / weight[g];
    for (unsigned int g = 0; g < n_groups; ++g)
      dst.m_sigf[g] = dst.m_sigf[g] / weight[g];

    // Renormalizing chi if necessary
    double total_chi = 0;
    for (unsigned int g = 0; g < n_groups; ++g)
      total_chi += dst.m_chi[g];
    const double eps = 1.e-16;
    if (std::abs (total_chi) > eps)
    {
      for (unsigned int g = 0; g < n_groups; ++g)
        dst.m_chi[g] = dst.m_chi[g] / total_chi;
    }
    else
    {
      dst.m_chi[0] = 1.0;
      for (unsigned int g = 1; g < n_groups; ++g)
        dst.m_chi[g] = 0;
    }

    dst.calculate_diff ();
    dst.calculate_sigmaa ();
    dst.calculate_sigmar ();
  }

  /*template <int dim>
   void
   Homogenization<dim>::print_homogenized_xs (const std::string &filename,
   const std::vector<
   InputMat::XS_single> & xs,
   const stdvd3 & phi0,
   const stdvd3 & phi2,
   const stdvd3 & current)
   {
   // Do they exist?
   const bool exist_phi0 = true;
   const bool exist_phi2 = true;
   const bool exist_current = true;

   // We need some constants
   const unsigned int n_groups = mp_state.get_n_groups ();

   // Create empty property tree object
   using boost::property_tree::ptree;
   ptree pt;

   std::string bin; // bin is a container to save temporal data
   std::string indent ("  ");

   // materials ---------------------------------------------------------
   pt.add ("materials", "");
   pt.put ("materials.<xmlattr>.ngroups", n_groups);

   for (unsigned int i = 0; i < xs.size (); ++i)
   {
   // here we generate the node
   ptree & node = pt.add ("materials.mix", "");

   node.put ("<xmlattr>.id", i);
   node.put ("name", "this is a pin");
   // node.put ("name", xs_mat->second.name);

   if (xs[i].exist_diff)
   {
   Forest::vector_to_string (xs[i].m_diff, bin);
   node.put ("Diff", bin);
   bin.clear ();
   }

   if (xs[i].exist_sigmat)
   {
   Forest::vector_to_string (xs[i].m_sigmat, bin);
   node.put ("SigmaT", bin);
   bin.clear ();
   }

   if (xs[i].exist_sigmaa)
   {
   Forest::vector_to_string (xs[i].m_sigmaa, bin);
   node.put ("SigmaA", bin);
   bin.clear ();
   }

   if (xs[i].exist_sigmar)
   {
   Forest::vector_to_string (xs[i].m_sigmar, bin);
   node.put ("SigmaR", bin);
   bin.clear ();
   }

   if (xs[i].exist_nusigf)
   {
   Forest::vector_to_string (xs[i].m_nusigf, bin);
   node.put ("NuSigF", bin);
   bin.clear ();
   }

   if (xs[i].exist_nu)
   {
   Forest::vector_to_string (xs[i].m_nu, bin);
   node.put ("Nu", bin);
   bin.clear ();
   }

   if (xs[i].exist_sigf)
   {
   Forest::vector_to_string (xs[i].m_sigf, bin);
   node.put ("SigF", bin);
   bin.clear ();
   }

   if (xs[i].exist_chi)
   {
   Forest::vector_to_string (xs[i].m_chi, bin);
   node.put ("Chi", bin);
   bin.clear ();
   }

   if (xs[i].exist_sigmas)
   {
   Forest::vector_to_string (xs[i].m_sigmas, bin, indent, 3);
   node.put ("SigmaS", bin);
   bin.clear ();
   }

   if (exist_phi0)
   {
   node.add ("<xmlcomment>", "Phi_0 per face and per group");
   Forest::vector_to_string (phi0[i], bin, indent, 3);
   node.put ("Phi_0", bin);
   bin.clear ();
   }

   if (exist_phi2)
   {
   node.add ("<xmlcomment>", "Phi_2 per face and per group");
   Forest::vector_to_string (phi2[i], bin, indent, 3);
   node.put ("Phi_2", bin);
   bin.clear ();
   }

   if (exist_current)
   {
   Forest::vector_to_string (current[i], bin, indent, 3);
   node.put ("Current", bin);
   bin.clear ();
   }

   }

   // Write property tree to XML file (\t is a tabulator)
   boost::property_tree::xml_writer_settings<std::string> settings (' ', 2);
   write_xml (filename, pt, std::locale (), settings);
   }*/

  template class Homogenization<1> ;
  template class Homogenization<2> ;
  template class Homogenization<3> ;

} // end of namespace Forest
