#ifndef _SECONDORDERLEAPFROG_HH_
#define _SECONDORDERLEAPFROG_HH_

#include "TimeIntegrator.hh"

// =================================================================================================
/** @brief The class SecondOrderLeapFrog.

    Second-order Kick-Drift-Kick (KDK) Leapfrog integration scheme.

    Within a single Move call, given acceleration a_n at time t_n:
      Half kick:  v_half    = v_n + (dt/2) * a_n
      Drift:      x_{n+1}  = x_n + dt * v_half         (=> x_n + dt*v_n + (dt^2/2)*a_n)
      Half kick:  v_{n+1}  = v_half + (dt/2) * a_n     (=> v_n + dt*a_n)

    The rotational analogue uses the same KDK structure:
      Half kick:  omega_half  = omega_n + (dt/2) * alpha_n
      Drift:      q_{n+1}    = q_n x dq(omega_half)
      Half kick:  omega_{n+1} = omega_half + (dt/2) * alpha_n

    The position/orientation update is second-order accurate in time. The
    velocity update is formally first-order when only a_n is available (the
    second half-kick using a_{n+1} would require the next force evaluation),
    but the drift step already provides the second-order spatial accuracy that
    distinguishes this scheme from first-order symplectic Euler.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
class SecondOrderLeapFrog : public TimeIntegrator<T>
{
public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    SecondOrderLeapFrog();

    /** @brief Constructor with the time step
        @param dt time step */
    __HOSTDEVICE__
    SecondOrderLeapFrog(T dt);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~SecondOrderLeapFrog();
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

    /** @brief Computes half-kick on velocity and full drift on position/orientation.
        This is Step 1 of KDK: v_{n+1/2} = v_n + (dt/2)*a_n, x_{n+1} = x_n + dt*v_{n+1/2}.
        After this call the stored velocity equals v_{n+1/2}.
        @param momentum acceleration at t_n
        @param velocity velocity (updated in-place to v_{n+1/2})
        @param transMotion translational displacement dt*v_{n+1/2}
        @param rotMotion rotational change (quaternion increment) */
    __HOSTDEVICE__
    void Move(const Kinematics<T>& momentum,
              Kinematics<T>&       velocity,
              Vector3<T>&          transMotion,
              Quaternion<T>&       rotMotion) const final;

    /** @brief Performs the second velocity half-kick: v_{n+1} = v_{n+1/2} + (dt/2)*a_{n+1}.
        This is Step 3 of KDK, called after forces at x_{n+1} have been computed.
        @param momentum acceleration at t_{n+1}
        @param velocity velocity (updated in-place from v_{n+1/2} to v_{n+1}) */
    __HOSTDEVICE__
    void AdvanceVelocity(const Kinematics<T>& momentum, Kinematics<T>& velocity) const final;
    //@}
};

#endif
