/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/hom_data.cc
 * @brief  Implementation of class template HomData
 */

#include "neutronics/hom_data.h"
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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION < 105500)
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#endif

#include <string>

namespace Forest
{
  using namespace dealii;

  template <int dim>
  HomData<dim>::HomData (const State<dim> & state,
                         const Homogenization<dim> & homogenization,
                         const ExtractBV<dim> & extract_bv)
      : mp_state (state),
        mp_homogenization (homogenization),
        mp_extract_bv (extract_bv),
        data (mp_state.mp_data),
        geom (mp_state.mp_geom)
  {
    // HomPinData--------------------------------------------------------------
    // name of the output file
    std::string f_pin_xs (
        data.mp_settings.get_path_plus_file () + ".pin.mat.xml");
    // We make this members of homogenization visible in this function.
    const std::vector<XS_single> & pin_xs = homogenization.pin_xs;
    // We make this members of extract_bv visible in this function.
    const stdvd4 & pin_moments = mp_extract_bv.m_pin_bvm;
    // now we print the values
    print_homogenized_xs (f_pin_xs, pin_xs, pin_moments);
    /*// We make this members of extract_bv visible in this function.
    const stdvd3 & phi0_pin = mp_extract_bv.m_phi0_pin;
    const stdvd3 & phi2_pin = mp_extract_bv.m_phi2_pin;
    const stdvd3 & current_pin = mp_extract_bv.m_j_pin;
    // now we print the values
    print_homogenized_xs (f_pin_xs, pin_xs, phi0_pin, phi2_pin, current_pin);*/

    // HomAssData--------------------------------------------------------------
    // name of the output file
    std::string f_ass_xs (
        data.mp_settings.get_path_plus_file () + ".ass.mat.xml");
    // We make this members of homogenization visible in this function.
    const std::vector<XS_single> & ass_xs = homogenization.ass_xs;
    // We make this members of extract_bv visible in this function.
    const stdvd4 & ass_moments = mp_extract_bv.m_ass_bvm;
    // now we print the values
    print_homogenized_xs (f_ass_xs, ass_xs, ass_moments);
    /*// We make this members of extract_bv visible in this function.
    const stdvd3 & phi0_ass = mp_extract_bv.m_phi0_ass;
    const stdvd3 & phi2_ass = mp_extract_bv.m_phi2_ass;
    const stdvd3 & current_ass = mp_extract_bv.m_j_ass;
    // now we print the values
    print_homogenized_xs (f_ass_xs, ass_xs, phi0_ass, phi2_ass, current_ass);*/
  }

  template <int dim>
  void
  HomData<dim>::print_homogenized_xs (const std::string &filename,
                                             const std::vector<
                                                 XS_single> & xs,
                                             const stdvd3 & phi0,
                                             const stdvd3 & phi2,
                                             const stdvd3 & current)
  {
    const bool exist_phi0 = true;
    const bool exist_phi2 = true;
    const bool exist_current = true;

    /* We need some constants */
    const unsigned int n_groups = mp_state.get_n_groups ();

    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    std::string bin; // bin is a container to save temporal data
    std::string indent ("  ");

    // materials ---------------------------------------------------------
    pt.add ("materials", "");
    pt.put ("materials.<xmlattr>.ngroups", n_groups);
    pt.add ("materials.<xmlcomment>", "K-effective for the reactor");
    pt.put ("materials.keff", mp_state.get_keff ());
    for (unsigned int i = 0; i < xs.size (); ++i)
    {
      // here we generate the node
      ptree & node = pt.add ("materials.mix", "");

      node.put ("<xmlattr>.id", i);
      node.put ("<xmlattr>.name", "Material");
      /* node.put ("name", xs_mat->second.name); */

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
    write_xml (filename, pt, std::locale (), xml_writer_settings(' ', 2));
  }

  template <int dim>
  void
  HomData<dim>::print_homogenized_xs (const std::string &filename,
                                             const std::vector<XS_single> & xs,
                                             const stdvd4 & phi_moments)
  {
    /* We need some constants */
    const unsigned int n_groups = mp_state.get_n_groups ();

    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    std::string bin; // bin is a container to save temporal data
    std::string indent ("  ");

    // materials ---------------------------------------------------------
    pt.add ("materials", "");
    pt.put ("materials.<xmlattr>.ngroups", n_groups);
    pt.add ("materials.<xmlcomment>", "K-effective for the reactor");
    pt.put ("materials.keff", mp_state.get_keff ());

    for (unsigned int i = 0; i < xs.size (); ++i)
    {
      // here we generate the node
      ptree & node = pt.add ("materials.mix", "");

      node.put ("<xmlattr>.id", i);
      node.put ("<xmlattr>.name", "Material");
      /* node.put ("name", xs_mat->second.name); */

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

      for (unsigned int m = 0; m < phi_moments[i].size(); ++m)
      {
        node.add ("<xmlcomment>", "Phi_"+num_to_string (m)+" per face and per group");
        Forest::vector_to_string (phi_moments[i][m], bin, indent, 3);
        node.put ("Phi_"+num_to_string (m), bin);
        bin.clear ();
      }

    }

    // Write property tree to XML file (\t is a tabulator)
    write_xml (filename, pt, std::locale (), xml_writer_settings(' ', 2));
  }

  template class HomData<1> ;
  template class HomData<2> ;
  template class HomData<3> ;

} // end of namespace Forest
