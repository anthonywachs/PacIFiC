#ifndef _GRAINSTESTDEV_HH_
#define _GRAINSTESTDEV_HH_

#include "Grains.hh"

// =================================================================================================
/** @brief The class GrainsTestDev.

    For developers to test implementations.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class GrainsTestDev : public Grains<T>
{
public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    GrainsTestDev();

    /** @brief Destructor */
    ~GrainsTestDev();
    //@}

    /** @name High-level methods */
    //@{
    /** @brief Runs the simulation over the prescribed time interval */
    void simulate() final;
    //@}
};

#endif
