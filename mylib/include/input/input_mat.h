/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   input/input_mat.h
 * @brief  InputMat class template declarations
 */
#ifndef FOREST_INPUT_MAT_H
#define FOREST_INPUT_MAT_H

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Forest
{

  /**
   * @brief Structure containing the data of a single material.
   * @details Different data defining the properties of a single material,
   * and booleans variable to check which data has been provided
   * to this material or not.
   */
  class XS_single
  {
  public:

    XS_single ()
        : n_groups (0),
          id (-1),
          exist_sigmat (false),
          exist_sigmas (false),
          exist_nusigf (false),
          exist_chi (false)
    {
    }
    ;

    XS_single (const unsigned int id,
               const std::vector<double> & st,
               const std::vector<std::vector<double> > & ss,
               const std::vector<double> & nsf,
               const std::vector<double> & chi,
               const std::string name = "xs")
        : n_groups (st.size ()),
          id (id),
          exist_sigmat (true),
          exist_sigmas (true),
          exist_nusigf (true),
          exist_chi (true)
    {
      m_sigmat = st;
      m_sigmas = ss;
      m_nusigf = nsf;
      m_chi = chi;
      this->name = name;
    }
    ;

    unsigned int n_groups;

    /** @brief Unique identification number for this material */
    unsigned int id;

    /** @brief Name of the material (just information). */
    std::string name;

    /** @brief flag saying that we have filled the cross section */
    bool exist_sigmat = false;
    /** @brief Total cross section for each energy group. */
    std::vector<double> m_sigmat;

    /** @brief flag saying that we have filled the cross section */
    bool exist_sigmas = false;
    /** @brief Scattering cross section matrix for each energy pair (g,h). */
    std::vector<std::vector<double> > m_sigmas;

    /** @brief flag saying that we have filled the cross section */
    bool exist_nusigf = false;
    /** @brief Fission cross section for each energy group. */
    std::vector<double> m_nusigf;

    /** @brief flag saying that we have filled the cross section */
    bool exist_chi = false;
    /** @brief Fission spectrum for each energy group. */
    std::vector<double> m_chi;

    /** @brief flag saying that we have filled the cross section */
    bool exist_sigmaa = false;
    /** @brief Absortion cross section for each energy group. */
    std::vector<double> m_sigmaa;

    /** @brief flag saying that we have filled the cross section */
    bool exist_nu = false;
    /** @brief fission cross section for each energy group. */
    std::vector<double> m_nu;

    /** @brief flag saying that we have filled the cross section */
    bool exist_sigf = false;
    /** @brief fission cross section for each energy group. */
    std::vector<double> m_sigf;

    /** @brief flag saying that we have filled the cross section */
    bool exist_diff = false;
    /** @brief Diffusion coefficient for each energy group. */
    std::vector<double> m_diff;

    /** @brief flag saying that we have filled the cross section */
    bool exist_sigmar = false;
    /** @brief Absortion cross section for each energy group. */
    std::vector<double> m_sigmar;

    /** @brief flag saying that we have filled the cross section */
    bool exist_pitch = false;
    /** @brief Absortion cross section for each energy group. */
    std::vector<double> m_pitch;

  public:
    /**
     * @brief Checks that `nusigf[g] == nu[g]*sigf[g]` with certain tolerance.
     * @details It assumes that all the cross sections exists.
     */
    void
    check_nusigf () const;

    /**
     * @brief Checks that `sigmat[g] == sigmaa[g] + sigmas[g][:]`.
     * @details It assumes that all the cross sections exists.
     */
    void
    check_sigmat () const;

    /** @brief Calculates `nusigf[g]` from `nu[g]*sigf[g]` */
    void
    calculate_nusigf ();

    /** @brief Calculates `sigmat[g]` from `sigmaa[g] + sigmas[g][:]`. */
    void
    calculate_sigmat ();

    /** @brief Calculates `sigmaa[g]` from `sigmat[g] - sigmas[g][:]`. */
    void
    calculate_sigmaa ();

    /** @brief Calculates `diff[g]` from `1/(.3*sigmat[g])`. */
    void
    calculate_diff ();

    /** @brief Calculates `sigmar[g]` from `sigmaa[g] + sigmas[g][1:g-1,g+1:]`. */
    void
    calculate_sigmar ();

    /** @brief normalize chi to sum equal to 1.  */
    void
    normalize_chi ();
  };

  /**
   * @class InputMat
   * @brief Class for the cross sections.
   * @ingroup ForestInput
   * @details Here we have the cross sections for all the materials,
   * as well as some functions to read and write this data
   * to xml format, and some functions to check the data
   * when enough cross sections are available, i.e.,
   * we check that
   * \f$ \nu\Sigma_{f} = \nu*\Sigma_{f} \f$
   *  and we check
   * \f$ \Sigma_{t,g} = \sum_{h}(\Sigma_{s,g,h}) +
   * \Sigma_{a,g} \f$
   */
  class InputMat
  {
  public:

    /** @brief Constructor */
    InputMat ();
    InputMat (const std::string &filename);

    /** @brief Short for the pair used for the cross sections */
    typedef std::pair<unsigned int, XS_single> XS_pair;
    /** @brief Short for the map used for the cross sections */
    typedef std::map<unsigned int, XS_single> XS_map;

    /** @brief Map to store the cross section data. */
    XS_map m_xs;

    void
    set (const XS_single & xs_mat);

    void
    set (const unsigned int id,
         const std::vector<double> & st,
         const std::vector<std::vector<double> > & ss,
         const std::vector<double> & nsf,
         const std::vector<double> & chi,
         const std::string name = "Mixture");

    /** @brief Load the material data from the file @p filename */
    void
    load (const std::string &filename);

    /**
     @brief Print the material data to the file @p filename
     */
    void
    save (const std::string &filename);

    /**
     @brief Check Consistency of the materials data.
     @details check the data
     when enough cross sections are available, i.e.,
     we check that
     \f$ \nu\Sigma_{f} = \nu*\Sigma_{f} \f$
     and we check
     \f$ \Sigma_{t,g} = \sum_{h}(\Sigma_{s,g,h}) + \Sigma_{a,g} \f$
     */
    void
    check ();

    /**
     * @brief Get the number of energy groups.
     * @return n_groups
     */
    unsigned int
    get_n_groups () const
    {
      return n_groups;
    }

    /**
     * @brief Get the number of different materials defined.
     * @return n_mat
     */
    unsigned int
    get_n_mat () const
    {
      return n_mat;
    }

    /**
     * @brief Return a pointer to the xs structure of the particular material
     */
    /*
     XS_single & const get_xs(unsigned int i_mat_) const
     {
     return *xs[i_mat_];
     }
     */

  private:

    /** Number of energy groups for the materials. */
    unsigned int n_groups;

    /** Number different materials. */
    unsigned int n_mat;

  };

} // end of namespace Forest

#endif /* FOREST_INPUT_MAT_H */
