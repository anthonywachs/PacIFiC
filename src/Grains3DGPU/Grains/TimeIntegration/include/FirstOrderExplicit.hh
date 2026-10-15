#ifndef _FIRSTORDEREXPLICIT_HH_
#define _FIRSTORDEREXPLICIT_HH_

#include "TimeIntegrator.hh"

// =================================================================================================
/** @brief The class FirstOrderExplicit.

    First order explicit integration scheme: x(t+dt) = x(t) + dt*v(t) and
    v(t+dt) = v(t) + dt*a(t).

    @author A.Yazdani - 2024 - Construction */
// ============================================================================
template <typename T>
class FirstOrderExplicit : public TimeIntegrator<T>
{
public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    FirstOrderExplicit();

    /** @brief Constructor with the time step */
    __HOSTDEVICE__
    FirstOrderExplicit(T dt);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~FirstOrderExplicit();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the time integrator type */
    __HOSTDEVICE__
    TimeIntegratorType getTimeIntegratorType() const final;
    //@}

    /** @name Methods */
    //@{
    /** @brief Creates and returns a clone of the time integrator */
    __HOSTDEVICE__
    TimeIntegrator<T>* clone() const final;

    /** @brief Computes the new velocity and transformation change over dt
        @param momentum acceleration
        @param velocity velocity
        @param transMotion translational motion over dt
        @param rotMotion rotational motion over dt */
    __HOSTDEVICE__
    void Move(const Kinematics<T>& momentum,
              Kinematics<T>&       velocity,
              Vector3<T>&          transMotion,
              Quaternion<T>&       rotMotion) const final;
    //@}
};

#endif
