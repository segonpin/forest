/**
 * @author Sebastian Gonzalez-Pintor. Chalmers University, 2014.
 * @file   neutronics/scattering_source.h
 * @brief  ScatteringSource class template declarations
 */

#ifndef FOREST_SCATTERING_SOURCE_H
#define FOREST_SCATTERING_SOURCE_H

#include "algebra/sn_vector_fwd.h"
#include "deal.II/lac/block_vector.h"   // for BlockVector

namespace Forest
{
  using namespace dealii;

  /**
   * @brief Scattering Source, base class
   * @ingroup ForestNeutronics
   */
  template <int dim>
  class ScatteringSource
  {

  public:

    /** @brief Virtual Destructor */
    virtual
    ~ScatteringSource ()
    {
    }

    /** @brief Return true when operator is matrix free, false otherwise. */
    virtual bool
    is_matrix_free() = 0;

  public:

    /** @todo do I need it?
    void
    vmult_add (SnVector & source,
           const SnVector & phi)
    {
      for (unsigned int g = 0; g < source.m_data.size (); ++g)
        for (unsigned int h = 0; h < phi.m_data.size (); ++h)
          vmult_add_g_h (source.m_data[g], phi.m_data[h], g, h);
    }

    void
    vmult (SnVector & source,
           const SnVector & phi)
    {
      source = 0;
      vmult_add (source, phi);
    }
    */

    /** @todo Document me */
    void
    down_scatt_add (SnVector & source,
                    const SnVector & phi,
                    const unsigned int g);

    /**
     * @todo Document me
     * @warning Here we assume that the size of the vector phi is the same
     * as the number of groups. It is true for the moment, but I should take
     * care about possible modifications of the snvector and don't keep the
     * implementation based on the assumption.
     */
    void
    up_scatt_add (SnVector & source,
                  const SnVector & phi,
                  const unsigned int g);

    /** @todo Document me */
    void
    in_scatt (SnVector & source,
              const SnVector & phi,
              const unsigned int g);

    /** @todo Document me */
    void
    in_scatt (BlockVector<double> & source,
              const BlockVector<double> & phi,
              const unsigned int g);

    /** @todo Document me */
    void
    in_scatt_add (SnVector & source,
                  const SnVector & phi,
                  const unsigned int g);

    /**
     @brief Calculates the memory consumption of this class.
     @return The memory consumption of this class. Zero for matrix free mode.
     */
    virtual double
    memory_consumption () const = 0;

    void
    vmult_add_g_h (SnVector & source,
                   const SnVector & phi,
                   const unsigned int g,
                   const unsigned int h);

    /**
     @brief This is the building block for the functions in the class.
     @details Applied in a group wise level, switch between matrix-free or
     non-matrix-free to calculate the scattering from g to h of phi and add it
     to the value of source.
     @param source
     @param phi
     @param g
     @param h
     */
    virtual void
    vmult_add_g_h (BlockVector<double> & source,
                   const BlockVector<double> & phi,
                   const unsigned int g,
                   const unsigned int h) = 0;


  public:
    virtual bool
    no_upscatt(const unsigned int g) const = 0;

  };

} // end of namespace Forest

#endif /* FOREST_SCATTERING_SOURCE_H */
