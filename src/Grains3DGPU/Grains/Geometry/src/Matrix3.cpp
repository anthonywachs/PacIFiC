#include "Matrix3.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor. Matrix is initialized to the identity matrix
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3() noexcept
{
    setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
}

// -------------------------------------------------------------------------------------------------
// Constructor with a 1D array of values as input
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3(T const* buffer) noexcept
{
    setValue(buffer);
}

// -------------------------------------------------------------------------------------------------
// Constructor with 9 components as inputs
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) noexcept
{
    setValue(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

// -------------------------------------------------------------------------------------------------
// Constructor with 3 angles (in radians) as input parameters)
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3(T roll, T pitch, T yaw) noexcept
{
    T cr = cos(roll);
    T sr = sin(roll);
    T cp = cos(pitch);
    T sp = sin(pitch);
    T cy = cos(yaw);
    T sy = sin(yaw);

    T xx = cy * cp;
    T xy = cy * sp * sr - sy * cr;
    T xz = cy * sp * cr + sy * sr;

    T yx = sy * cp;
    T yy = sy * sp * sr + cy * cr;
    T yz = sy * sp * cr - cy * sr;

    T zx = -sp;
    T zy = cp * sr;
    T zz = cp * cr;

    setValue(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3(const Matrix3<T>& mat) noexcept
{
    setValue(mat.getBuffer());
}

// -------------------------------------------------------------------------------------------------
// Assign operator to another matrix
template <typename T>
__HOSTDEVICE__ Matrix3<T>& Matrix3<T>::operator=(const Matrix3<T>& m) noexcept
{
    if(&m != this)
        setValue(m.getBuffer());
    return (*this);
}

// -------------------------------------------------------------------------------------------------
// Move constructor
template <typename T>
__HOSTDEVICE__ Matrix3<T>::Matrix3(Matrix3<T>&& mat) noexcept
{
    setValue(mat.getBuffer());
    mat.setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
}

// -------------------------------------------------------------------------------------------------
// Move assignment operator
template <typename T>
__HOSTDEVICE__ Matrix3<T>& Matrix3<T>::operator=(Matrix3<T>&& m) noexcept
{
    if(&m != this)
    {
        setValue(m.getBuffer());
        m.setValue(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
    }
    return (*this);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ Matrix3<T>::Matrix3(DOMNode* root) noexcept
{
    if(root)
    {
        std::string        values = ReaderXML::getNodeValue_String(root);
        std::istringstream inValues(values.c_str());
        inValues >> this->m_comp[XX] >> this->m_comp[XY] >> this->m_comp[XZ] >> this->m_comp[YX]
            >> this->m_comp[YY] >> this->m_comp[YZ] >> this->m_comp[ZX] >> this->m_comp[ZY]
            >> this->m_comp[ZZ];
    }
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Matrix3<T>::~Matrix3() noexcept
{
}

// -------------------------------------------------------------------------------------------------
/* Gets the pointer to the buffer */
template <typename T>
__HOSTDEVICE__ T const* Matrix3<T>::getBuffer() const noexcept
{
    return (m_comp);
}

// -------------------------------------------------------------------------------------------------
// Sets the matrix to a 1D array of 9 values as input
template <typename T>
__HOSTDEVICE__ void Matrix3<T>::setValue(T const* buffer) noexcept
{
    m_comp[XX] = buffer[XX];
    m_comp[XY] = buffer[XY];
    m_comp[XZ] = buffer[XZ];
    m_comp[YX] = buffer[YX];
    m_comp[YY] = buffer[YY];
    m_comp[YZ] = buffer[YZ];
    m_comp[ZX] = buffer[ZX];
    m_comp[ZY] = buffer[ZY];
    m_comp[ZZ] = buffer[ZZ];
}

// -------------------------------------------------------------------------------------------------
// Sets the matrix with all 9 components as inputs
template <typename T>
__HOSTDEVICE__ void
    Matrix3<T>::setValue(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) noexcept
{
    m_comp[XX] = xx;
    m_comp[XY] = xy;
    m_comp[XZ] = xz;
    m_comp[YX] = yx;
    m_comp[YY] = yy;
    m_comp[YZ] = yz;
    m_comp[ZX] = zx;
    m_comp[ZY] = zy;
    m_comp[ZZ] = zz;
}

// -------------------------------------------------------------------------------------------------
// i-th row accessor
template <typename T>
__HOSTDEVICE__ Vector3<T>& Matrix3<T>::operator[](uint i) noexcept
{
    return (*(Vector3<T>*)(m_comp + 3 * i));
}

// -------------------------------------------------------------------------------------------------
// i-th row accessor
template <typename T>
__HOSTDEVICE__ const Vector3<T>& Matrix3<T>::operator[](uint i) const noexcept
{
    return (*(const Vector3<T>*)(m_comp + 3 * i));
}

// -------------------------------------------------------------------------------------------------
// element accessor
template <typename T>
__HOSTDEVICE__ T& Matrix3<T>::operator()(uint i) noexcept
{
    return (m_comp[i]);
}

// -------------------------------------------------------------------------------------------------
// const element accessor
template <typename T>
__HOSTDEVICE__ const T& Matrix3<T>::operator()(uint i) const noexcept
{
    return (m_comp[i]);
}

// -------------------------------------------------------------------------------------------------
// element accessor
template <typename T>
__HOSTDEVICE__ T& Matrix3<T>::operator()(uint i, uint j) noexcept
{
    return (m_comp[i * 3 + j]);
}

// -------------------------------------------------------------------------------------------------
// const element accessor
template <typename T>
__HOSTDEVICE__ const T& Matrix3<T>::operator()(uint i, uint j) const noexcept
{
    return (m_comp[i * 3 + j]);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Matrix3<T>& m)
{
    fileOut << m[X] << "\n" << m[Y] << "\n" << m[Z];
    return (fileOut);
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Matrix3<T>& m)
{
    fileIn >> m[X] >> m[Y] >> m[Z];
    return (fileIn);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Matrix3<float>;
template class Matrix3<double>;

#define X(T)                                                                            \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Matrix3<T>& m); \
                                                                                        \
    template std::istream& operator>> <T>(std::istream & fileIn, Matrix3<T> & m);
X(float)
X(double)
#undef X
