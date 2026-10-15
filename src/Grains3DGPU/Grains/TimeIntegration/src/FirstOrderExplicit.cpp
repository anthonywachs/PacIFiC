#include "FirstOrderExplicit.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ FirstOrderExplicit<T>::FirstOrderExplicit()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with the time step
template <typename T>
__HOSTDEVICE__ FirstOrderExplicit<T>::FirstOrderExplicit(T dt)
{
    TimeIntegrator<T>::m_dt = dt;
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ FirstOrderExplicit<T>::~FirstOrderExplicit()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the time integrator type
template <typename T>
__HOSTDEVICE__ TimeIntegratorType FirstOrderExplicit<T>::getTimeIntegratorType() const
{
    return (FIRSTORDEREXPLICIT);
}

// -------------------------------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
template <typename T>
__HOSTDEVICE__ TimeIntegrator<T>* FirstOrderExplicit<T>::clone() const
{
    return (new FirstOrderExplicit<T>(TimeIntegrator<T>::m_dt));
}

// -------------------------------------------------------------------------------------------------
// Computes the new velocity and transformation change over dt
template <typename T>
__HOSTDEVICE__ void FirstOrderExplicit<T>::Move(const Kinematics<T>& momentum,
                                                Kinematics<T>&       velocity,
                                                Vector3<T>&          transMotion,
                                                Quaternion<T>&       rotMotion) const
{
    T dt = TimeIntegrator<T>::m_dt;
    // Semi-implicit (symplectic) Euler for better stability with stiff forces:
    //   v_{n+1} = v_n + dt * a_n
    //   x_{n+1} = x_n + dt * v_{n+1}
    velocity.addToTranslationalComponent(dt * momentum.getTranslationalComponent());
    transMotion = dt * velocity.getTranslationalComponent();

    velocity.addToAngularComponent(dt * momentum.getAngularComponent());
    rotMotion = this->computeQuaternionChange(velocity.getAngularComponent());
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class FirstOrderExplicit<float>;
template class FirstOrderExplicit<double>;
