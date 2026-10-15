#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ Torce<T>::Torce()
    : m_torque(0, 0, 0)
    , m_force(0, 0, 0)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with a torque and a force as input parameters
template <typename T>
__HOSTDEVICE__ Torce<T>::Torce(const Vector3<T>& t, const Vector3<T>& f)
    : m_torque(t)
    , m_force(f)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Torce<T>::~Torce()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the total torque of the torce
template <typename T>
__HOSTDEVICE__ Vector3<T> Torce<T>::getTorque() const
{
    return (m_torque);
}

// -------------------------------------------------------------------------------------------------
// Gets the total force of the torce
template <typename T>
__HOSTDEVICE__ Vector3<T> Torce<T>::getForce() const
{
    return (m_force);
}

// -------------------------------------------------------------------------------------------------
// Sets the total torque of the torce
template <typename T>
__HOSTDEVICE__ void Torce<T>::setTorque(const Vector3<T>& t)
{
    m_torque = t;
}

// -------------------------------------------------------------------------------------------------
// Sets the total force of the torce
template <typename T>
__HOSTDEVICE__ void Torce<T>::setForce(const Vector3<T>& f)
{
    m_force = f;
}

// -------------------------------------------------------------------------------------------------
// Resets the torce
template <typename T>
__HOSTDEVICE__ void Torce<T>::reset()
{
    m_torque.reset();
    m_force.reset();
}

// -------------------------------------------------------------------------------------------------
// Adds a force to the torce
template <typename T>
__HOSTDEVICE__ void Torce<T>::addTorque(const Vector3<T>& t)
{
    m_torque += t;
}

// -------------------------------------------------------------------------------------------------
// Adds a force to the torce
template <typename T>
__HOSTDEVICE__ void Torce<T>::addForce(const Vector3<T>& f)
{
    m_force += f;
}

// -------------------------------------------------------------------------------------------------
// Adds a force to the torce with accounting for the additional torque
template <typename T>
__HOSTDEVICE__ void Torce<T>::addForce(const Vector3<T>& f, const Vector3<T>& p)
{
    m_force += f;
    m_torque += (p ^ f);
}

// -------------------------------------------------------------------------------------------------
// Atomically adds a torque to the torce (GPU safe)
template <typename T>
__device__ void Torce<T>::addTorqueAtomic(const Vector3<T>& t)
{
    m_torque.atomicAdd(t);
}

// -------------------------------------------------------------------------------------------------
// Atomically adds a force to the torce (GPU safe)
template <typename T>
__device__ void Torce<T>::addForceAtomic(const Vector3<T>& f)
{
    m_force.atomicAdd(f);
}

// -------------------------------------------------------------------------------------------------
// Atomically adds both force and torque from another Torce (GPU safe)
template <typename T>
__device__ void Torce<T>::addTorceAtomic(const Torce<T>& other)
{
    m_force.atomicAdd(other.m_force);
    m_torque.atomicAdd(other.m_torque);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Torce<T>& t)
{
    fileOut << "Torque: " << t.getTorque() << "\n"
            << "Force: " << t.getForce();
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Torce<T>& t)
{
    Vector3<T> vec;
    fileIn >> vec;
    t.setTorque(vec);
    fileIn >> vec;
    t.setForce(vec);
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Torce<float>;
template class Torce<double>;

#define X(T)                                                                          \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Torce<T>& t); \
    template std::istream& operator>> <T>(std::istream & fileIn, Torce<T> & t);
X(float)
X(double)
#undef X