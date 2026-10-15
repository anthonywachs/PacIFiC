#include "SecondOrderLeapFrog.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ SecondOrderLeapFrog<T>::SecondOrderLeapFrog()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with the time step
template <typename T>
__HOSTDEVICE__ SecondOrderLeapFrog<T>::SecondOrderLeapFrog(T dt)
{
    TimeIntegrator<T>::m_dt = dt;
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ SecondOrderLeapFrog<T>::~SecondOrderLeapFrog()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the time integrator type
template <typename T>
__HOSTDEVICE__ TimeIntegratorType SecondOrderLeapFrog<T>::getTimeIntegratorType() const
{
    return (SECONDORDERLEAPFROG);
}

// -------------------------------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
template <typename T>
__HOSTDEVICE__ TimeIntegrator<T>* SecondOrderLeapFrog<T>::clone() const
{
    return (new SecondOrderLeapFrog<T>(TimeIntegrator<T>::m_dt));
}

// -------------------------------------------------------------------------------------------------
// Step 1 of KDK: half-kick velocity then drift position/orientation
template <typename T>
__HOSTDEVICE__ void SecondOrderLeapFrog<T>::Move(const Kinematics<T>& momentum,
                                                 Kinematics<T>&       velocity,
                                                 Vector3<T>&          transMotion,
                                                 Quaternion<T>&       rotMotion) const
{
    const T dt      = TimeIntegrator<T>::m_dt;
    const T half_dt = T(0.5) * dt;

    // -- Translational --
    // Half kick: v_half = v_n + (dt/2) * a_n
    const Vector3<T> vHalf
        = velocity.getTranslationalComponent() + half_dt * momentum.getTranslationalComponent();
    // Drift: dx = dt * v_half  (second-order position update)
    transMotion = dt * vHalf;
    // Store v_{n+1/2} -- second half-kick happens later in AdvanceVelocity
    velocity.setTranslationalComponent(vHalf);

    // -- Rotational --
    // Half kick: omega_half = omega_n + (dt/2) * alpha_n
    const Vector3<T> omegaHalf
        = velocity.getAngularComponent() + half_dt * momentum.getAngularComponent();
    // Drift: dq computed from omega_half (second-order orientation update)
    rotMotion = this->computeQuaternionChange(omegaHalf);
    // Store omega_{n+1/2}
    velocity.setAngularComponent(omegaHalf);
}

// -------------------------------------------------------------------------------------------------
// Step 3 of KDK: second half-kick using forces at the new position
template <typename T>
__HOSTDEVICE__ void SecondOrderLeapFrog<T>::AdvanceVelocity(const Kinematics<T>& momentum,
                                                            Kinematics<T>&       velocity) const
{
    const T half_dt = T(0.5) * TimeIntegrator<T>::m_dt;
    // v_{n+1} = v_{n+1/2} + (dt/2) * a_{n+1}
    velocity.addToTranslationalComponent(half_dt * momentum.getTranslationalComponent());
    velocity.addToAngularComponent(half_dt * momentum.getAngularComponent());
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class SecondOrderLeapFrog<float>;
template class SecondOrderLeapFrog<double>;
