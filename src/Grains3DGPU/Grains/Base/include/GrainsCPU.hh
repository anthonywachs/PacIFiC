#ifndef _GRAINSCPU_HH_
#define _GRAINSCPU_HH_

#include "Grains.hh"

// =================================================================================================
/** @brief The class GrainsCPU.

    Standard Grains3D application on CPU.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class GrainsCPU : public Grains<T>
{
protected:
    /** @name Parameters */
    //@{
    //@}

public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    GrainsCPU();

    /** @brief Destructor */
    ~GrainsCPU();
    //@}

    /** @name High-level methods */
    //@{
    // /** @brief Tasks to perform before time-stepping
    // @param rootElement XML root */
    // void initialize( DOMElement* rootElement );

    /** @brief Runs the simulation over the prescribed time interval */
    void simulate() final;

    // /** @brief Tasks to perform after time-stepping */
    // virtual void finalize();
    //@}
};

#endif
