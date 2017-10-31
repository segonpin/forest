/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input_mat.cc
 * @brief  Implementation of class InputMat
 */

#include "input/input_mat.h"
#include "utils/constants.h"

#include "utils/forest_utils_base.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION < 105500)
  typedef boost::property_tree::xml_writer_settings<char> xml_writer_settings;
#else
  typedef boost::property_tree::xml_writer_settings<std::string> xml_writer_settings;
#endif
//#include <boost/foreach.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>                      // for basic_filebuf<>::int_type, etc
#include <sstream>                      // for basic_stringbuf<>::int_type, etc

#include <deal.II/base/exceptions.h>    // for Assert and ExcMessage
#include <cmath>

namespace Forest
{

  InputMat::InputMat ()
      : n_groups (0),
        n_mat (0)
  {
  }

  InputMat::InputMat (const std::string &filename)
      : n_groups (0),
        n_mat (0)
  {
    // We read the cross sections from the input file
    this->load (filename);

    // Now we check that the cross sections are consistent, and generate
    // the missing ones
    this->check ();

    //this->save(input_path + problem_name+"_out.mat.xml");
  }

  void
  InputMat::set (const XS_single & xs_mat)
  {
    m_xs.insert (XS_pair (xs_mat.id, xs_mat));
    //std::cout << "set material" << xs_mat.id << '\n';
    // update n_groups if different
    if (n_groups != xs_mat.n_groups)
    {
      // if different from zero, it means it was already initialized
      // but it has changed
      // (error??, different materials with different number of groups??)
      if (n_groups != 0)
        std::cout << "warning: changing the number of energy groups with " <<
        xs_mat.name << " and mat.id = " << xs_mat.id << '\n';
      n_groups = xs_mat.n_groups;
    }
    // then update the number of different materials
    n_mat = m_xs.size ();
    //std::cout << "Now we have " << n_mat << " materials and " << n_groups << " energy groups" << '\n';
  }

  void
  InputMat::set (const unsigned int id,
                 const std::vector<double> & st,
                 const std::vector<std::vector<double> > & ss,
                 const std::vector<double> & nsf,
                 const std::vector<double> & chi,
                 const std::string name)
  {
    //XS_single xs_mat(id, st, ss, nsf, chi, name);
    //this->set(xs_mat);
    this->set(XS_single(id, st, ss, nsf, chi, name));
  }


  void
  InputMat::load (const std::string &filename)
  {
    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    // Load XML file and put its contents in property tree. 
    // No namespace qualification is needed, because of Koenig 
    // lookup on the second argument. If reading fails, exception
    // is thrown.
    read_xml (filename, pt, boost::property_tree::xml_parser::trim_whitespace);

    n_groups = pt.get<unsigned int> ("materials.<xmlattr>.ngroups");

    // Range based loop (c++11) and auto (c++11) must be the way to go.
    // for (const auto &mat : pt.get_child("materials")
    // for (const std::pair<std::string, ptree> & mat : pt.get_child("materials"))
    // { if (mat.first != "mix") continue; }

    // cross sections ----------------------------------------------------
    for (const std::pair<std::string, ptree> & v : pt.get_child ("materials"))
    {
      // If it is not a mixture, then we jump it
      if (v.first != "mix")
        continue;

      XS_single xs_mat;
      xs_mat.n_groups = n_groups;
      xs_mat.id = v.second.get<unsigned int> ("<xmlattr>.id");
      xs_mat.name = v.second.get<std::string> ("<xmlattr>.name",
          std::string ("mix"));

      std::string bin; // bin is a container to store temporal data

      bin = v.second.get<std::string> ("SigmaT", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_sigmat);
        bin.clear ();
        xs_mat.exist_sigmat = true;
      }

      bin = v.second.get<std::string> ("SigmaS", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_sigmas);
        bin.clear ();
        xs_mat.exist_sigmas = true;
      }

      bin = v.second.get<std::string> ("SigmaA", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_sigmaa);
        bin.clear ();
        xs_mat.exist_sigmaa = true;
      }

      bin = v.second.get<std::string> ("NuSigF", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_nusigf);
        bin.clear ();
        xs_mat.exist_nusigf = true;
      }

      bin = v.second.get<std::string> ("Nu", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_nu);
        bin.clear ();
        xs_mat.exist_nu = true;
      }

      bin = v.second.get<std::string> ("SigF", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_sigf);
        bin.clear ();
        xs_mat.exist_sigf = true;
      }

      bin = v.second.get<std::string> ("Chi", std::string (""));
      if (bin != std::string (""))
      {
        Forest::string_to_vector (bin, xs_mat.m_chi);
        bin.clear ();
        xs_mat.exist_chi = true;
      }

      m_xs.insert (XS_pair (xs_mat.id, xs_mat));
    }
    n_mat = m_xs.size ();
  }

  void
  InputMat::save (const std::string &filename)
  {
    // Create empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    std::string bin; // bin is a container to save temporal data
    std::string indent ("\t");

    // materials ---------------------------------------------------------
    pt.add ("materials", "");
    pt.put ("materials.<xmlattr>.ngroups", n_groups);

    for (typename XS_map::iterator xs = m_xs.begin (); xs != m_xs.end (); ++xs)
    {
      // here we generate the node
      ptree & node = pt.add ("materials.mix", "");

      node.put ("<xmlattr>.id", xs->first);
      node.put ("<xmlattr>.name", xs->first);

      if (xs->second.exist_sigmat)
      {
        Forest::vector_to_string (xs->second.m_sigmat, bin);
        node.put ("SigmaT", bin);
        bin.clear ();
      }

      if (xs->second.exist_sigmaa)
      {
        Forest::vector_to_string (xs->second.m_sigmaa, bin);
        node.put ("SigmaA", bin);
        bin.clear ();
      }

      if (xs->second.exist_nusigf)
      {
        Forest::vector_to_string (xs->second.m_nusigf, bin);
        node.put ("NuSigF", bin);
        bin.clear ();
      }

      if (xs->second.exist_nu)
      {
        Forest::vector_to_string (xs->second.m_nu, bin);
        node.put ("Nu", bin);
        bin.clear ();
      }

      if (xs->second.exist_sigf)
      {
        Forest::vector_to_string (xs->second.m_sigf, bin);
        node.put ("SigF", bin);
        bin.clear ();
      }

      if (xs->second.exist_chi)
      {
        Forest::vector_to_string (xs->second.m_chi, bin);
        node.put ("Chi", bin);
        bin.clear ();
      }

      if (xs->second.exist_sigmas)
      {
        Forest::vector_to_string (xs->second.m_sigmas, bin, indent, 3);
        node.put ("SigmaS", bin);
        bin.clear ();
      }

    }

    // Write property tree to XML file (\t is a tabulator)
    write_xml (filename, pt, std::locale (), xml_writer_settings (' ', 1));
  }

  void
  XS_single::check_nusigf () const
  {
    const double nusigf_reltol = 1.e-3;
    const double nusigf_abstol = 1.e-5;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double nusigf = m_nu[g] * m_sigf[g];
      if (std::abs (nusigf) < nusigf_reltol)
      {
        if (std::abs (nusigf - m_nusigf[g]) > nusigf_abstol)
        {
          std::cout << "Material id = " << id
                    << " has a different nusigf cross-section than the product"
                    << " of its nu and fission cross-sections for group " << g
                    << " : nusigf = " << std::setprecision (6) << m_nusigf[g]
                    << " calc_nusigf = " << std::setprecision (6) << nusigf
                    << "\n";
          //Assert(false,
          //    dealii::ExcMessage("Cross sections fail consistency check."));
        }
      }
      else
      {
        if (std::abs ((nusigf - m_nusigf[g]) / nusigf) > nusigf_reltol)
        {
          std::cout << "Material id = " << id
                    << " has a different nusigf cross-section than the product"
                    << " of its nu and fission cross-sections for group " << g
                    << " : nusigf = " << std::setprecision (6) << m_nusigf[g]
                    << " calc_nusigf = " << std::setprecision (6) << nusigf
                    << "\n";
          //Assert(false,
          //    dealii::ExcMessage("Cross section fail consistency check."));
        }
      }
    }
  }

  void
  XS_single::check_sigmat () const
  {
    const double sigmat_threshold = 1.e-4;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double sigmat = m_sigmaa[g];
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        sigmat += m_sigmas[h][g];
      }
      if (std::abs (sigmat - m_sigmat[g]) > sigmat_threshold)
      {
        std::cout
            << "Material id = " << id
            << " has a different total cross-section than the sum"
            << " of its scattering and absorption cross-sections for group "
            << g << " : sigma_t = " << std::setprecision (6) << m_sigmat[g]
            << " calc_sigma_t = " << std::setprecision (6) << sigmat << "\n";
        //Assert(false,
        //    dealii::ExcMessage("Cross sections fail consistency check."));
      }
    }
  }

  void
  XS_single::calculate_nusigf ()
  {
    Assert(exist_nu and exist_sigf,
        dealii::ExcMessage("Prerequisites not satisfied"));
    m_nusigf.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double nusigf = m_nu[g] * m_sigf[g];
      m_nusigf[g] = nusigf;
    }
    exist_nusigf = true;
  }

  void
  XS_single::calculate_sigmat ()
  {
    Assert(exist_sigmaa and exist_sigmas,
        dealii::ExcMessage("Prerequisites not satisfied"));
    m_sigmat.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double sigmat = m_sigmaa[g];
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        sigmat += m_sigmas[h][g];
      }
      m_sigmat[g] = sigmat;
    }
    exist_sigmat = true;
  }

  void
  XS_single::calculate_sigmaa ()
  {
    Assert(exist_sigmat and exist_sigmas,
        dealii::ExcMessage("Prerequisites not satisfied"));
    m_sigmaa.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double sigmaa = m_sigmat[g];
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        sigmaa -= m_sigmas[h][g];
      }
      m_sigmaa[g] = sigmaa;
    }
    exist_sigmaa = true;
  }

  void
  XS_single::calculate_diff ()
  {
    Assert(exist_sigmat, dealii::ExcMessage("Prerequisites not satisfied"));
    m_diff.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      Assert(std::abs(m_sigmat[g])>1.e-6,
          dealii::ExcMessage("Zero total xs generates infinite diffusion."));
      m_diff[g] = 1. / (3. * m_sigmat[g]);
    }
    exist_diff = true;
  }

  void
  XS_single::calculate_sigmar ()
  {
    Assert(exist_sigmat and exist_sigmas,
        dealii::ExcMessage("Prerequisites not satisfied"));
    m_sigmar.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      m_sigmar[g] = m_sigmat[g] - m_sigmas[g][g];
    }
    exist_sigmar = true;
    /*Assert(exist_sigmaa and exist_sigmas,
        dealii::ExcMessage("Prerequisites not satisfied"));
    m_sigmar.resize (n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double sigmar = m_sigmaa[g] - m_sigmas[g][g];
      for (unsigned int h = 0; h < n_groups; ++h)
      {
        sigmar += m_sigmas[h][g];
      }
      m_sigmar[g] = sigmar;
    }
    exist_sigmar = true;*/
  }

  void
  XS_single::normalize_chi ()
  {
    Assert(exist_chi,
        dealii::ExcMessage("Prerequisites not satisfied"));
    double sum_chi = 0.;
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      sum_chi += m_chi[g];
    }
    Assert(std::abs(sum_chi) > eps6,
        dealii::ExcMessage("Wrong chi adds to zero."));
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      m_chi[g] = m_chi[g] / sum_chi;
    }
  }

  void
  InputMat::check ()
  {
    for (typename XS_map::iterator xs = m_xs.begin (); xs != m_xs.end (); ++xs)
    {
      // Transport and Diffusion cross sections -------------------------------

      // Check for nusigf and calculate it in case it is not provided.
      if (xs->second.exist_nu && !xs->second.exist_sigf)
      {
        // in case it was provided, check consistency for given tolerance.
        if (xs->second.exist_nusigf)
        {
          xs->second.check_nusigf ();
        }
        xs->second.calculate_nusigf ();
      }

      // Normalize chi in case it was provided (should always be provided?)
      if (xs->second.exist_chi)
      {
        xs->second.normalize_chi ();
      }

      // Transport cross sections ---------------------------------------------

      // Check consistency for sigmat in case several data is provided
      if (xs->second.exist_sigmas && xs->second.exist_sigmaa)
      {
        // Check consistency for sigmat in case it is provided
        if (xs->second.exist_sigmat)
        {
          xs->second.check_sigmat ();
        }
        // Calculate sigmat, provided or not.
        xs->second.calculate_sigmat ();
      }
      else if (xs->second.exist_sigmas && xs->second.exist_sigmat)
      {
        // sigmas and sigmat exist, but not sigmaa, so we calculate it
        xs->second.calculate_sigmaa ();
      }

      // Diffusion cross sections ---------------------------------------------

      if (xs->second.exist_sigmat && xs->second.exist_sigmas)
        xs->second.calculate_diff ();

      if (xs->second.exist_sigmaa && xs->second.exist_sigmas)
        xs->second.calculate_sigmar ();

#ifdef DEBUG

      // Check if the problem has been completely defined for transport
      bool problem_defined = false;
      if (xs->second.exist_sigmat && xs->second.exist_sigmas
          && xs->second.exist_chi && xs->second.exist_nusigf)
      {
        problem_defined = true;
      }
      // Check if the problem has been completely defined for diffusion
      if (xs->second.exist_diff && xs->second.exist_sigmar
          && xs->second.exist_sigmas && xs->second.exist_chi
          && xs->second.exist_nusigf)
      {
        problem_defined = true;
      }
      Assert(problem_defined, dealii::ExcMessage("Problem not fully defined."));
#endif
    }

  }

} // end of namespace Forest
