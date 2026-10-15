#include "TimeIntegrator.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ TimeIntegrator<T>::TimeIntegrator()
{
}

// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ TimeIntegrator<T>::TimeIntegrator(TimeIntegrator<T> const& ti)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ TimeIntegrator<T>::~TimeIntegrator()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the time step
template <typename T>
__HOSTDEVICE__ T TimeIntegrator<T>::getTimeStep() const
{
    return m_dt;
}

// -------------------------------------------------------------------------------------------------
// Computes the quaternion change over the time step
template <typename T>
__HOSTDEVICE__ Quaternion<T>
               TimeIntegrator<T>::computeQuaternionChange(const Vector3<T>& avgAngVel) const
{
    // Quaternion change over dt
    const T nOmega = norm(avgAngVel);
    if(nOmega > EPS<T>)
    {
        const T c = cos(nOmega * m_dt / T(2));
        const T s = sin(nOmega * m_dt / T(2));
        // T c = cos( nOmega * m_dt );
        // T s = sin( nOmega * m_dt );
        return (Quaternion<T>((s / nOmega) * avgAngVel, c));
    }
    else
        return (Quaternion<T>(T(0), T(1)));
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class TimeIntegrator<float>;
template class TimeIntegrator<double>;