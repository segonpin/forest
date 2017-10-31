#ifndef REAL_CORE_H
#define REAL_CORE_H

#include "input/input_fwd.h"            // For class Input&
#include "input/input_geom.h"    // for Core

#include "deal.II/base/point.h"  // for Point

#include <vector>                // for vector

#ifndef DOXYGEN
namespace Forest { template <int dim> class RealAssembly; }
namespace Forest { template <int dim> class RealPin; }
namespace dealii { template <int dim, int spacedim> class Triangulation; }
#endif

namespace Forest
{
  using namespace dealii;

  /// this name is provisional. Something better can be better.
  template <int dim>
  class RealCore
  {
  public:

    unsigned int m_empty_assembly = static_cast<unsigned int> (-1);

    /**
     * @brief Constructor using the input data
     * @param input_data_
     * @param correction_factor
     */
    RealCore (const Input & input_data_,
              const double correction_factor = 1.0);

    /**
     *
     */
    void
    build_data ();

    /**
     *
     * @param tria
     */
    void
    build_mesh (Triangulation<dim, dim> & tria);

    /**
     * @todo Study if we needed. If not remove it
     * @param correction_factor
     */
    inline void
    set_correction_factor (const unsigned int correction_factor)
    {
      m_correction_factor = correction_factor;
    }

    /**
     *
     * @return
     */
    inline unsigned int
    get_n_assemblies () const
    {
      return m_n_assemblies;
    }

    /**
     *
     * @return
     */
    inline unsigned int
    get_n_pins () const
    {
      return m_n_pins;
    }

    /**
     *
     * @return
     */
    inline std::vector<unsigned int>
    get_pins_to_assembly () const
    {
      return pins_to_assembly;
    }

    /**
     *
     * @return
     */
    inline std::vector<unsigned int>
    get_pins_per_assembly () const
    {
      return pins_per_assembly;
    }

    /**
     *
     * @param ass_no
     * @return
     */
    inline unsigned int
    get_n_pins_in_assembly (const unsigned int ass_no) const
    {
      return pins_per_assembly[ass_no];
    }

  private:
    /**
     *
     */
    void
    build_assemblies ();

    /**
     *
     */
    void
    build_additional_data ();

    /// @todo I should change this to private
  public:
    const Input & m_input_data;

    typename Forest::Core m_core;

    Point<dim> m_origin;

    /**
     * @brief  size is a vector with the size
     *  (sum?? of all the pin sizes) of the assembly
     */
    std::vector<double> m_size;

    std::vector<RealAssembly<dim> > m_real_assemblies;
    std::vector<RealPin<dim> > m_real_pins;

    bool m_ready_to_build_mesh;
    double m_correction_factor;

    unsigned int m_n_assemblies;
    unsigned int m_n_pins;

    std::vector<unsigned int> pins_to_assembly;
    std::vector<unsigned int> pins_per_assembly;
  };
}

#endif

