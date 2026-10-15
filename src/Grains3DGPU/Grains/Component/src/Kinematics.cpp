#include "Kinematics.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ Kinematics<T>::Kinematics()
    : m_translational(0, 0, 0)
    , m_angular(0, 0, 0)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with two vectors as input parameters
template <typename T>
__HOSTDEVICE__ Kinematics<T>::Kinematics(const Vector3<T>& translational, const Vector3<T>& angular)
    : m_translational(translational)
    , m_angular(angular)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Kinematics<T>::~Kinematics()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the translational component of the kinematics
template <typename T>
__HOSTDEVICE__ Vector3<T> Kinematics<T>::getTranslationalComponent() const
{
    return (m_translational);
}

// -------------------------------------------------------------------------------------------------
// Gets the angular component of the kinematics
template <typename T>
__HOSTDEVICE__ Vector3<T> Kinematics<T>::getAngularComponent() const
{
    return (m_angular);
}

// -------------------------------------------------------------------------------------------------
// Sets the translational component of the kinematics
template <typename T>
__HOSTDEVICE__ void Kinematics<T>::setTranslationalComponent(const Vector3<T>& translational)
{
    m_translational = translational;
}

// -------------------------------------------------------------------------------------------------
// Sets the angular component of the kinematics
template <typename T>
__HOSTDEVICE__ void Kinematics<T>::setAngularComponent(const Vector3<T>& angular)
{
    m_angular = angular;
}

// -------------------------------------------------------------------------------------------------
// Adds a vector to the translational component of the kinematics
template <typename T>
__HOSTDEVICE__ void Kinematics<T>::addToTranslationalComponent(const Vector3<T>& translational)
{
    m_translational += translational;
}

// -------------------------------------------------------------------------------------------------
// Adds a vector to the angular component of the kinematics
template <typename T>
__HOSTDEVICE__ void Kinematics<T>::addToAngularComponent(const Vector3<T>& omega)
{
    m_angular += omega;
}

// -------------------------------------------------------------------------------------------------
// Returns the total velocity U + om x R given R
template <typename T>
__HOSTDEVICE__ Vector3<T> Kinematics<T>::kinematicsAtPoint(const Vector3<T>& R) const
{
    return (m_translational + (m_angular ^ R));
    // return ( m_translational );
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Kinematics<T>& k)
{
    fileOut << "Translational Component: " << k.getTranslationalComponent() << "\n"
            << "Angular Component: " << k.getAngularComponent();
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Kinematics<T>& k)
{
    Vector3<T> vec;
    fileIn >> vec;
    k.setTranslationalComponent(vec);
    fileIn >> vec;
    k.setAngularComponent(vec);
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Kinematics<float>;
template class Kinematics<double>;

#define X(T)                                                                               \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Kinematics<T>& k); \
                                                                                           \
    template std::istream& operator>> <T>(std::istream & fileIn, Kinematics<T> & k);
X(float)
X(double)
#undef X