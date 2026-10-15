#ifndef _TIMEINTEGRATORFACTORY_HH_
#define _TIMEINTEGRATORFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "TimeIntegrator.hh"

// =================================================================================================
/** @brief The class TimeIntegratorFactory.

    Creates the numerical scheme for the time integration of the Newton's law
    and the kinematic equations.

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring
    @author A.YAZDANI - 2024 - Major cleaning for porting to GPU */
// =================================================================================================
template <typename T>
class TimeIntegratorFactory
{
private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    TimeIntegratorFactory();

    /** @brief Destructor (forbidden) */
    ~TimeIntegratorFactory();
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates and returns the time integration scheme.
        @param root XML node
        @param dt time step
        @param TI time integrator buffer */
    __HOST__
    static void create(DOMNode* root, T dt, GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>& TI);

    /** @brief TimeIntegrator objects must be instantiated on device, if we want to use them on
        device. Copying from host is not supported due to runtime polymorphism for this class. This
        function reads a host-side TimeIntegrator object, and mimics it in a given device buffer.
        It calls a device kernel that is implemented in the source file.
        @param h_TI Host-side TimeIntegrator object
        @param d_TI Device-side TimeIntegrator object */
    __HOST__
    static void copyHostToDevice(GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>&   h_TI,
                                 GrainsMemBuffer<TimeIntegrator<T>*, MemType::DEVICE>& d_TI);
    //@}
};

#endif
