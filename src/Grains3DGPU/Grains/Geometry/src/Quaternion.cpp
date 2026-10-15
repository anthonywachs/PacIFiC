#include "Quaternion.hh"
#include "GrainsUtils.hh"
#include "MatrixMath.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion() noexcept
    : m_vqt(T(0))
    , m_w(T(1))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with 2 scalar as input parameters q and d. Quaternion is
// initialized as [ d, (q,q,q) ]
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(T q, T w) noexcept
    : m_vqt(q)
    , m_w(w)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with a Vector3 vector vec and a scalar d. Quaternion is
// initialized as [ d, vec ]
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(const Vector3<T>& vec, T w) noexcept
    : m_vqt(vec)
    , m_w(w)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with a vector given by its 3 components (x,y,z) and a scalar d.
// Quaternion is initialized as [ d, (x,y,z) ]
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(T x, T y, T z, T w) noexcept
    : m_vqt(Vector3<T>(x, y, z))
    , m_w(w)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with a buffer
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(const T* buffer) noexcept
    : m_vqt(Vector3<T>(buffer))
    , m_w(buffer[3])
{
}

// -------------------------------------------------------------------------------------------------
// Constructor from Euler angles (Z-Y-X intrinsic order)
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(T aX, T aY, T aZ) noexcept
{
    setQuaternion(aX, aY, aZ);
}

// -------------------------------------------------------------------------------------------------
// Constructor with a rotation matrix
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(const Matrix3<T>& rot) noexcept
{
    this->setQuaternion(rot);
}

// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(const Quaternion<T>& q) noexcept
    : m_vqt(q.m_vqt)
    , m_w(q.m_w)
{
}

// -------------------------------------------------------------------------------------------------
// Assign operator to another Quaternion object
template <typename T>
__HOSTDEVICE__ Quaternion<T>& Quaternion<T>::operator=(const Quaternion<T>& q) noexcept
{
    if(this != &q)  // self-assignment check
    {
        m_vqt = q.m_vqt;
        m_w   = q.m_w;
    }
    return *this;
}

// -------------------------------------------------------------------------------------------------
// Move constructor
template <typename T>
__HOSTDEVICE__ Quaternion<T>::Quaternion(Quaternion<T>&& q) noexcept
    : m_vqt(std::move(q.m_vqt))
    , m_w(q.m_w)
{
    q.m_w = T(0);
    q.m_vqt.reset();
}

// -------------------------------------------------------------------------------------------------
// Move assignment operator
template <typename T>
__HOSTDEVICE__ Quaternion<T>& Quaternion<T>::operator=(Quaternion<T>&& q) noexcept
{
    if(this != &q)  // self-assignment check
    {
        m_vqt = std::move(q.m_vqt);
        m_w   = q.m_w;
        q.m_w = T(0);
        q.m_vqt.reset();
    }
    return *this;
}

// -------------------------------------------------------------------------------------------------
// Constructor from an XML node
template <typename T>
__HOST__ Quaternion<T>::Quaternion(DOMNode* root) noexcept
{
    std::string type = ReaderXML::getNodeAttr_String(root, "Type");
    if(type == "Matrix")
    {
        Matrix3<T> mat(root);
        setQuaternion(mat);
    }
    else if(type == "Angles")
    {
        // read in radiands
        T aX = RADS_PER_DEG<T> * T(ReaderXML::getNodeAttr_Double(root, "aX"));
        T aY = RADS_PER_DEG<T> * T(ReaderXML::getNodeAttr_Double(root, "aY"));
        T aZ = RADS_PER_DEG<T> * T(ReaderXML::getNodeAttr_Double(root, "aZ"));

        setQuaternion(aX, aY, aZ);
    }
    else if(type == "Identity")
    {
        setQuaternion(Vector3<T>(0, 0, 0), T(1));
    }
    else
        GAbort("A quaternion in one of the AngularPosition XML nodes is"
               " not a rotation matrix or angle.");
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Quaternion<T>::~Quaternion() noexcept
{
}

// -------------------------------------------------------------------------------------------------
// Returns the pointer to the buffer
// The buffer is an array of 4 elements: [x, y, z, w]
template <typename T>
__HOSTDEVICE__ const T* Quaternion<T>::getBuffer() const noexcept
{
    return (reinterpret_cast<const T*>(&m_vqt));
}

// -------------------------------------------------------------------------------------------------
// Returns the vectorial part of the quaternion
template <typename T>
__HOSTDEVICE__ const Vector3<T>& Quaternion<T>::getVector() const noexcept
{
    return (m_vqt);
}

// -------------------------------------------------------------------------------------------------
// Returns the value of the scalar part of the quaternion
template <typename T>
__HOSTDEVICE__ const T& Quaternion<T>::getScalar() const noexcept
{
    return (m_w);
}

// -------------------------------------------------------------------------------------------------
// Sets the vectorial part of the quaternion
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setVector(const Vector3<T>& vec) noexcept
{
    m_vqt = vec;
}

// -------------------------------------------------------------------------------------------------
// Sets the scalar part of the quaternion
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setScalar(const T w) noexcept
{
    m_w = w;
}

// -------------------------------------------------------------------------------------------------
// Sets the quaternion with a Vector3 vector vec and a scalar d.
// Quaternion is set to [ d, vec ]
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setQuaternion(const Vector3<T>& vec, const T w) noexcept
{
    m_vqt = vec;
    m_w   = w;
}

// -------------------------------------------------------------------------------------------------
// Sets the quaternion with a vector given by its 3 components (x,y,z)
// and a scalar d. Quaternion is set to [ d, (x,y,z) ]
template <typename T>
__HOSTDEVICE__ void
    Quaternion<T>::setQuaternion(const T x, const T y, const T z, const T w) noexcept
{
    m_vqt[X] = x;
    m_vqt[Y] = y;
    m_vqt[Z] = z;
    m_w      = w;
}

// -------------------------------------------------------------------------------------------------
// Sets the quaternion with a rotation matrix
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setQuaternion(const Matrix3<T>& rot) noexcept
{
    // Validate that input is a proper rotation matrix
    GAssert(isRotation(rot),
            "Input matrix is not a valid rotation matrix in "
            "Quaternion::setQuaternion!");

    const T* b   = rot.getBuffer();  // rotation matrix buffer
    T        den = T(0);

    // Case rotYY > - rotZZ, rotXX > - rotYY and rotXX > - rotZZ
    if(b[YY] > -b[ZZ] && b[XX] > -b[YY] && b[XX] > -b[ZZ])
    {
        den = sqrt(T(1) + b[XX] + b[YY] + b[ZZ]);
        GAssert(den > EPS<T>,
                "Numerical instability in Quaternion::setQuaternion - "
                "denominator too small!");
        m_w      = T(0.5) * den;
        m_vqt[X] = T(0.5) * (b[ZY] - b[YZ]) / den;
        m_vqt[Y] = T(0.5) * (b[XZ] - b[ZX]) / den;
        m_vqt[Z] = T(0.5) * (b[YX] - b[XY]) / den;
    }
    // Case rotYY < - rotZZ, rotXX > rotYY and rotXX > rotZZ
    else if(b[YY] < -b[ZZ] && b[XX] > b[YY] && b[XX] > b[ZZ])
    {
        den = sqrt(T(1) + b[XX] - b[YY] - b[ZZ]);
        GAssert(den > EPS<T>,
                "Numerical instability in Quaternion::setQuaternion - "
                "denominator too small!");
        m_w      = T(0.5) * (b[ZY] - b[YZ]) / den;
        m_vqt[X] = T(0.5) * den;
        m_vqt[Y] = T(0.5) * (b[XY] + b[YX]) / den;
        m_vqt[Z] = T(0.5) * (b[ZX] + b[XZ]) / den;
    }
    // Case rotYY > rotZZ, rotXX < rotYY and rotXX < - rotZZ
    else if(b[YY] > b[ZZ] && b[XX] < b[YY] && b[XX] < -b[ZZ])
    {
        den = sqrt(T(1) - b[XX] + b[YY] - b[ZZ]);
        GAssert(den > EPS<T>,
                "Numerical instability in Quaternion::setQuaternion - "
                "denominator too small!");
        m_w      = T(0.5) * (b[XZ] - b[ZX]) / den;
        m_vqt[X] = T(0.5) * (b[XY] + b[YX]) / den;
        m_vqt[Y] = T(0.5) * den;
        m_vqt[Z] = T(0.5) * (b[YZ] + b[ZY]) / den;
    }
    // Case rotYY < rotZZ, rotXX < - rotYY and rotXX < rotZZ
    else if(b[YY] < b[ZZ] && b[XX] < -b[YY] && b[XX] < b[ZZ])
    {
        den = sqrt(T(1) - b[XX] - b[YY] + b[ZZ]);
        GAssert(den > EPS<T>,
                "Numerical instability in Quaternion::setQuaternion - "
                "denominator too small!");
        m_w      = T(0.5) * (b[YX] - b[XY]) / den;
        m_vqt[X] = T(0.5) * (b[ZX] + b[XZ]) / den;
        m_vqt[Y] = T(0.5) * (b[YZ] + b[ZY]) / den;
        m_vqt[Z] = T(0.5) * den;
    }
    else
        GAbort("Case not covered in Quaternion::setQuaternion!");
}

// -------------------------------------------------------------------------------------------------
// Sets the quaternion from three angles in radians
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setQuaternion(T aX, T aY, T aZ) noexcept
{
    // Build rotation matrix
    Matrix3<T> mat(aX, aY, aZ);
    setQuaternion(mat);
}

// -------------------------------------------------------------------------------------------------
// Build a unit quaternion representing the rotation from u to v.
// The input vectors need not be normalised. */
// TODO: if the input vectors aren't normalized, normalize them and warn the
// user.
template <typename T>
__HOSTDEVICE__ void Quaternion<T>::setRotFromTwoVectors(const Vector3<T>& u,
                                                        const Vector3<T>& v) noexcept
{
    T          norm_u_norm_v = sqrt((u * u) * (v * v));
    T          real_part     = norm_u_norm_v + u * v;
    Vector3<T> vect;

    if(real_part < 1.e-6 * norm_u_norm_v)
    {
        /* If u and v are exactly opposite, rotate 180 degrees
        around an arbitrary orthogonal axis. Axis normalisation
        can happen later, when we normalise the quaternion. */
        real_part = T(0);
        vect      = fabs(u[0]) > fabs(u[2]) ? Vector3<T>(-u[1], u[0], T(0))
                                            : Vector3<T>(T(0), -u[2], u[1]);
    }
    else
    {
        /* Otherwise, build quaternion the standard way. */
        vect = u ^ v;
    }

    Quaternion<T> qq(vect[0], vect[1], vect[2], real_part);
    *this = (T(1) / norm(qq)) * qq;
}

// -------------------------------------------------------------------------------------------------
// Builds a matrix from the quaternion
template <typename T>
__HOSTDEVICE__ Matrix3<T> Quaternion<T>::toMatrix() const noexcept
{
    T x2 = m_vqt[X] + m_vqt[X];
    T y2 = m_vqt[Y] + m_vqt[Y];
    T z2 = m_vqt[Z] + m_vqt[Z];
    T xx = m_vqt[X] * x2;
    T xy = m_vqt[X] * y2;
    T xz = m_vqt[X] * z2;
    T yy = m_vqt[Y] * y2;
    T yz = m_vqt[Y] * z2;
    T zz = m_vqt[Z] * z2;
    T wx = m_w * x2;
    T wy = m_w * y2;
    T wz = m_w * z2;
    return (Matrix3<T>(T(1) - (yy + zz),
                       xy - wz,
                       xz + wy,
                       xy + wz,
                       T(1) - (xx + zz),
                       yz - wx,
                       xz - wy,
                       yz + wx,
                       T(1) - (xx + yy)));
}

// -------------------------------------------------------------------------------------------------
// Multiplies the quaternion on the right by another quaternion rhs,
// i.e., perform this x rhs, and return the vectorial part of this x rhs
template <typename T>
__HOSTDEVICE__ Vector3<T> Quaternion<T>::multToVector3(const Quaternion<T>& q) const noexcept
{
    Vector3<T> vtmp((m_vqt ^ q.m_vqt) + (m_w * q.m_vqt) + (q.m_w * m_vqt));
    return (vtmp);
}

// -------------------------------------------------------------------------------------------------
// ith-component accessor: (0,1,2) for the vector components and 3 forthe scalar
template <typename T>
__HOSTDEVICE__ T Quaternion<T>::operator[](size_t i) const noexcept
{
    return (i == 3 ? m_w : m_vqt[i]);
}

// -------------------------------------------------------------------------------------------------
// ith-component accessor: (0,1,2) for the vector components and 3 for the
// scalar - modifiable lvalue
template <typename T>
__HOSTDEVICE__ T& Quaternion<T>::operator[](size_t i) noexcept
{
    return (i == 3 ? m_w : m_vqt[i]);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
std::ostream& operator<<(std::ostream& fileOut, const Quaternion<T>& q)
{
    fileOut << q.getScalar() << " " << q.getVector();
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
std::istream& operator>>(std::istream& fileIn, Quaternion<T>& q)
{
    Vector3<T> vec;
    T          scalar;
    fileIn >> scalar >> vec;
    q.setScalar(scalar);
    q.setVector(vec);
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Quaternion<float>;
template class Quaternion<double>;

#define X(T)                                                                               \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Quaternion<T>& q); \
                                                                                           \
    template std::istream& operator>> <T>(std::istream & fileIn, Quaternion<T> & q);
X(float)
X(double)
#undef X