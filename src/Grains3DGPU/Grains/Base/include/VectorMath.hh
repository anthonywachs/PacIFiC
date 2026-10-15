#ifndef _VECTORMATH_HH_
#define _VECTORMATH_HH_

#include "Vector3.hh"

// =================================================================================================
/** @brief Miscellaneous Vector3 functions and operators as header-only.

    Defining important vector functions and operators here as static functions.
    It will increase the binary size, but the performance gain is much more
    appreciated.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name Vector3 math functions and operators */
//@{
/** @brief Returns the norm of the vector
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE T norm(const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    return (sqrt(buffer[0] * buffer[0] + buffer[1] * buffer[1] + buffer[2] * buffer[2]));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the norm squared of the vector
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE T norm2(const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    return (buffer[0] * buffer[0] + buffer[1] * buffer[1] + buffer[2] * buffer[2]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Determines if the vector is approximately zero or not
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE bool isApproxZero(const Vector3<T>& v, T tol = EPS<T>) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    return (fabs(buffer[0]) < tol && fabs(buffer[1]) < tol && fabs(buffer[2]) < tol);
}

// -------------------------------------------------------------------------------------------------
/** @brief Rounds the components of the vector to +-tol
    @param v the vector
    @param tol tolerance -- EPS defined in Basic.hh is the default */
template <typename T>
__HOSTDEVICE__ static INLINE void round(Vector3<T>& v, T tol = EPS<T>) noexcept
{
    T* __RESTRICT__ buffer = const_cast<T*>(v.getBuffer());
    for(uint i = 0; i < 3; ++i)
    {
        if(fabs(buffer[i]) < tol)
            buffer[i] = T(0);
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Vectors addition
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator+(const Vector3<T>& v1,
                                                  const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1 = v1.getBuffer();
    const T* __RESTRICT__ b2 = v2.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b1[i] + b2[i];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Vectors addition in-place
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T>& operator+=(Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(v1.getBuffer());
    const T* __RESTRICT__ b2 = v2.getBuffer();
    for(uint i = 0; i < 3; ++i)
        b1[i] += b2[i];
    return v1;
}

// -------------------------------------------------------------------------------------------------
/** @brief Vectors subtraction
@param v1 1st vector
@param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator-(const Vector3<T>& v1,
                                                  const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1 = v1.getBuffer();
    const T* __RESTRICT__ b2 = v2.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b1[i] - b2[i];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Vectors subtraction in-place
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T>& operator-=(Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(v1.getBuffer());
    const T* __RESTRICT__ b2 = v2.getBuffer();
    for(uint i = 0; i < 3; ++i)
        b1[i] -= b2[i];
    return v1;
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-vector multiplication
    @param d the multiplication factor
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator*(const T d, const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = d * buffer[i];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-vector multiplication in-place
    @param v the vector
    @param d the multiplication factor */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T>& operator*=(Vector3<T>& v, const T d) noexcept
{
    T* __RESTRICT__ buffer = const_cast<T*>(v.getBuffer());
    for(uint i = 0; i < 3; ++i)
        buffer[i] *= d;
    return v;
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar division
    @param d division factor
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator/(const Vector3<T>& v, const T d) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = buffer[i] / d;
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar division in-place
    @param v the vector
    @param d division factor */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T>& operator/=(Vector3<T>& v, const T d) noexcept
{
    T* __RESTRICT__ buffer = const_cast<T*>(v.getBuffer());
    for(uint i = 0; i < 3; ++i)
        buffer[i] /= d;
    return v;
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-vector dot product
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE T operator*(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1  = v1.getBuffer();
    const T* __RESTRICT__ b2  = v2.getBuffer();
    T                     out = T(0);
    for(uint i = 0; i < 3; ++i)
        out += b1[i] * b2[i];
    return (out);
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-vector cross product
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator^(const Vector3<T>& v1,
                                                  const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1 = v1.getBuffer();
    const T* __RESTRICT__ b2 = v2.getBuffer();
    T __RESTRICT__        out[3];
    out[0] = b1[1] * b2[2] - b1[2] * b2[1];
    out[1] = b1[2] * b2[0] - b1[0] * b2[2];
    out[2] = b1[0] * b2[1] - b1[1] * b2[0];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-vector cross product in-place
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T>& operator^=(Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(v1.getBuffer());
    const T* __RESTRICT__ b2 = v2.getBuffer();
    T __RESTRICT__        out[3];
    out[0] = b1[1] * b2[2] - b1[2] * b2[1];
    out[1] = b1[2] * b2[0] - b1[0] * b2[2];
    out[2] = b1[0] * b2[1] - b1[1] * b2[0];
    for(uint i = 0; i < 3; ++i)
        b1[i] = out[i];
    return v1;
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector equality operator
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator==(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1 = v1.getBuffer();
    const T* __RESTRICT__ b2 = v2.getBuffer();
    for(int i = 0; i < 3; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector inequality operator
    @param v1 1st vector
    @param v2 2nd vector */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator!=(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
    const T* __RESTRICT__ b1 = v1.getBuffer();
    const T* __RESTRICT__ b2 = v2.getBuffer();
    for(int i = 0; i < 3; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return true;
    }
    return false;
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector sign flip
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator-(const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ buffer = v.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = -buffer[i];
    return (Vector3<T>(out));
}
//@}

#endif