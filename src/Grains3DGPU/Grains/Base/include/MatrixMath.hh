#ifndef _MATRIXMATH_HH_
#define _MATRIXMATH_HH_

#include "Matrix3.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief Miscellaneous Matrix3 functions and operators as header-only.

    Defining important matrix functions and operators here as static functions.
    It will increase the binary size, but the performance gain is much more
    appreciated.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name Matrix3 math functions and operators */
//@{
/** @brief Matrix absolute
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> fabs(const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b = m.getBuffer();
    return (Matrix3<T>(fabs(b[XX]),
                       fabs(b[XY]),
                       fabs(b[XZ]),
                       fabs(b[YX]),
                       fabs(b[YY]),
                       fabs(b[YZ]),
                       fabs(b[ZX]),
                       fabs(b[ZY]),
                       fabs(b[ZZ])));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix absolute in-place
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void fabs(Matrix3<T>& m) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(m.getBuffer());
    b[XX]             = fabs(b[XX]);
    b[XY]             = fabs(b[XY]);
    b[XZ]             = fabs(b[XZ]);
    b[YX]             = fabs(b[YX]);
    b[YY]             = fabs(b[YY]);
    b[YZ]             = fabs(b[YZ]);
    b[ZX]             = fabs(b[ZX]);
    b[ZY]             = fabs(b[ZY]);
    b[ZZ]             = fabs(b[ZZ]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix determinant
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE T determinant(const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b    = m.getBuffer();
    T                     out0 = b[XX] * (b[YY] * b[ZZ] - b[YZ] * b[ZY]);
    T                     out1 = b[XY] * (b[YZ] * b[ZX] - b[YX] * b[ZZ]);
    T                     out2 = b[XZ] * (b[YX] * b[ZY] - b[YY] * b[ZX]);
    return (out0 + out1 + out2);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix transposition
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> transpose(const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b = m.getBuffer();
    return (Matrix3<T>(b[XX], b[YX], b[ZX], b[XY], b[YY], b[ZY], b[XZ], b[YZ], b[ZZ]));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix transposition in-place
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void transpose(Matrix3<T>& m) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(m.getBuffer());
    T               temp[3];
    temp[0] = b[XY];
    temp[1] = b[XZ];
    temp[2] = b[YZ];
    b[XY]   = b[YX];
    b[XZ]   = b[ZX];
    b[YZ]   = b[ZY];
    b[YX]   = temp[0];
    b[ZX]   = temp[1];
    b[ZY]   = temp[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix inverse
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> inverse(const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b = m.getBuffer();
    T __RESTRICT__        out[9];

    // Calculate cofactor matrix
    out[XX] = (b[YY] * b[ZZ] - b[YZ] * b[ZY]);
    out[XY] = (b[XZ] * b[ZY] - b[XY] * b[ZZ]);
    out[XZ] = (b[XY] * b[YZ] - b[XZ] * b[YY]);
    out[YX] = (b[YZ] * b[ZX] - b[YX] * b[ZZ]);
    out[YY] = (b[XX] * b[ZZ] - b[XZ] * b[ZX]);
    out[YZ] = (b[XZ] * b[YX] - b[XX] * b[YZ]);
    out[ZX] = (b[YX] * b[ZY] - b[YY] * b[ZX]);
    out[ZY] = (b[XY] * b[ZX] - b[XX] * b[ZY]);
    out[ZZ] = (b[XX] * b[YY] - b[XY] * b[YX]);

    // Calculate determinant
    T det = b[XX] * out[XX] + b[XY] * out[YX] + b[XZ] * out[ZX];
    if(fabs(det) < EPS<T>)
        printf("Matrix is not inversible!\n");

    // Scale by inverse determinant
    T s = T(1) / det;
    for(int i = 0; i < 9; ++i)
        out[i] *= s;

    return (Matrix3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix inverse in-place
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void inverse(Matrix3<T>& m) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(m.getBuffer());
    T __RESTRICT__  out[9];

    // Calculate cofactor matrix
    out[XX] = (b[YY] * b[ZZ] - b[YZ] * b[ZY]);
    out[XY] = (b[XZ] * b[ZY] - b[XY] * b[ZZ]);
    out[XZ] = (b[XY] * b[YZ] - b[XZ] * b[YY]);
    out[YX] = (b[YZ] * b[ZX] - b[YX] * b[ZZ]);
    out[YY] = (b[XX] * b[ZZ] - b[XZ] * b[ZX]);
    out[YZ] = (b[XZ] * b[YX] - b[XX] * b[YZ]);
    out[ZX] = (b[YX] * b[ZY] - b[YY] * b[ZX]);
    out[ZY] = (b[XY] * b[ZX] - b[XX] * b[ZY]);
    out[ZZ] = (b[XX] * b[YY] - b[XY] * b[YX]);

    // Calculate determinant
    T det = b[XX] * out[XX] + b[XY] * out[YX] + b[XZ] * out[ZX];
    if(fabs(det) < EPS<T>)
        printf("Matrix is not inversible!\n");

    // Scale by inverse determinant
    T s = T(1) / det;
    for(int i = 0; i < 9; ++i)
        out[i] *= s;

    m.setValue(out);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix scale
    @param m the matrix
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> scale(const Matrix3<T>& m, const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ b1 = m.getBuffer();
    const T* __RESTRICT__ b2 = v.getBuffer();
    return Matrix3<T>(b1[XX] * b2[0],
                      b1[XY] * b2[1],
                      b1[XZ] * b2[2],
                      b1[YX] * b2[0],
                      b1[YY] * b2[1],
                      b1[YZ] * b2[2],
                      b1[ZX] * b2[0],
                      b1[ZY] * b2[1],
                      b1[ZZ] * b2[2]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix scale in-place
    @param m the matrix
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE void scale(Matrix3<T>& m, const Vector3<T>& v) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(m.getBuffer());
    const T* __RESTRICT__ b2 = v.getBuffer();
    b1[XX] *= b2[0];
    b1[XY] *= b2[1];
    b1[XZ] *= b2[2];
    b1[YX] *= b2[0];
    b1[YY] *= b2[1];
    b1[YZ] *= b2[2];
    b1[ZX] *= b2[0];
    b1[ZY] *= b2[1];
    b1[ZZ] *= b2[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix addition
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> operator+(const Matrix3<T>& m1,
                                                  const Matrix3<T>& m2) noexcept
{
    const T* __RESTRICT__ b1 = m1.getBuffer();
    const T* __RESTRICT__ b2 = m2.getBuffer();
    T __RESTRICT__        out[9];
    for(uint i = 0; i < 9; ++i)
        out[i] = b1[i] + b2[i];
    return (Matrix3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix addition in-place
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void operator+=(Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(m1.getBuffer());
    const T* __RESTRICT__ b2 = m2.getBuffer();
    for(uint i = 0; i < 9; ++i)
        b1[i] += b2[i];
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix subtraction
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> operator-(const Matrix3<T>& m1,
                                                  const Matrix3<T>& m2) noexcept
{
    const T* __RESTRICT__ b1 = m1.getBuffer();
    const T* __RESTRICT__ b2 = m2.getBuffer();
    T __RESTRICT__        out[9];
    for(uint i = 0; i < 9; ++i)
        out[i] = b1[i] - b2[i];
    return (Matrix3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix subtraction in-place
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void operator-=(Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(m1.getBuffer());
    const T* __RESTRICT__ b2 = m2.getBuffer();
    for(uint i = 0; i < 9; ++i)
        b1[i] -= b2[i];
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-matrix product
    @param c the scalar
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> operator*(T c, const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b = m.getBuffer();
    T __RESTRICT__        out[9];
    for(uint i = 0; i < 9; ++i)
        out[i] = c * b[i];
    return (Matrix3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-matrix product in-place
    @param c the scalar
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Matrix3<T>& m, T c) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(m.getBuffer());
    for(uint i = 0; i < 9; ++i)
        b[i] *= c;
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix-vector product
    @param m the matrix
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator*(const Matrix3<T>& m, const Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ b1 = m.getBuffer();
    const T* __RESTRICT__ b2 = v.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b1[3 * i] * b2[0] + b1[3 * i + 1] * b2[1] + b1[3 * i + 2] * b2[2];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix-vector product in-place. Note that this modifies the vector.
    @param m the matrix
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(const Matrix3<T>& m, Vector3<T>& v) noexcept
{
    const T* __RESTRICT__ b1 = m.getBuffer();
    T* __RESTRICT__       b2 = const_cast<T*>(v.getBuffer());
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b1[3 * i] * b2[0] + b1[3 * i + 1] * b2[1] + b1[3 * i + 2] * b2[2];
    v.setValue(out);
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-matrix product
    @param v the vector
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator*(const Vector3<T>& v, const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b1 = v.getBuffer();
    const T* __RESTRICT__ b2 = m.getBuffer();
    T __RESTRICT__        out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b2[i] * b1[0] + b2[i + 3] * b1[1] + b2[i + 6] * b1[2];
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-matrix product in-place. Note that this modifies the vector.
    @param v the vector
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Vector3<T>& v, const Matrix3<T>& m) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(v.getBuffer());
    T const* __RESTRICT__ b2 = m.getBuffer();
    T                     out[3];
    for(uint i = 0; i < 3; ++i)
        out[i] = b1[0] * b2[i] + b1[1] * b2[i + 3] + b1[2] * b2[i + 6];
    v.setValue(out);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix-matrix product
@param m right matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> operator*(const Matrix3<T>& m1,
                                                  const Matrix3<T>& m2) noexcept
{
    T const* __RESTRICT__ b1 = m1.getBuffer();
    T const* __RESTRICT__ b2 = m2.getBuffer();
    return (Matrix3<T>(b1[XX] * b2[XX] + b1[XY] * b2[YX] + b1[XZ] * b2[ZX],
                       b1[XX] * b2[XY] + b1[XY] * b2[YY] + b1[XZ] * b2[ZY],
                       b1[XX] * b2[XZ] + b1[XY] * b2[YZ] + b1[XZ] * b2[ZZ],
                       b1[YX] * b2[XX] + b1[YY] * b2[YX] + b1[YZ] * b2[ZX],
                       b1[YX] * b2[XY] + b1[YY] * b2[YY] + b1[YZ] * b2[ZY],
                       b1[YX] * b2[XZ] + b1[YY] * b2[YZ] + b1[YZ] * b2[ZZ],
                       b1[ZX] * b2[XX] + b1[ZY] * b2[YX] + b1[ZZ] * b2[ZX],
                       b1[ZX] * b2[XY] + b1[ZY] * b2[YY] + b1[ZZ] * b2[ZY],
                       b1[ZX] * b2[XZ] + b1[ZY] * b2[YZ] + b1[ZZ] * b2[ZZ]));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix-matrix product in-place
    @param m1 left matrix
    @param m2 right matrix */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(m1.getBuffer());
    const T* __RESTRICT__ b2 = m2.getBuffer();
    T __RESTRICT__        out[9];
    out[XX] = b1[XX] * b2[XX] + b1[XY] * b2[YX] + b1[XZ] * b2[ZX];
    out[XY] = b1[XX] * b2[XY] + b1[XY] * b2[YY] + b1[XZ] * b2[ZY];
    out[XZ] = b1[XX] * b2[XZ] + b1[XY] * b2[YZ] + b1[XZ] * b2[ZZ];
    out[YX] = b1[YX] * b2[XX] + b1[YY] * b2[YX] + b1[YZ] * b2[ZX];
    out[YY] = b1[YX] * b2[XY] + b1[YY] * b2[YY] + b1[YZ] * b2[ZY];
    out[YZ] = b1[YX] * b2[XZ] + b1[YY] * b2[YZ] + b1[YZ] * b2[ZZ];
    out[ZX] = b1[ZX] * b2[XX] + b1[ZY] * b2[YX] + b1[ZZ] * b2[ZX];
    out[ZY] = b1[ZX] * b2[XY] + b1[ZY] * b2[YY] + b1[ZZ] * b2[ZY];
    out[ZZ] = b1[ZX] * b2[XZ] + b1[ZY] * b2[YZ] + b1[ZZ] * b2[ZZ];
    m1.setValue(out);
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix sign flip
    @param m the matrix */
template <typename T>
__HOSTDEVICE__ static INLINE Matrix3<T> operator-(const Matrix3<T>& m) noexcept
{
    const T* __RESTRICT__ b = m.getBuffer();
    return (Matrix3<T>(-b[XX], -b[XY], -b[XZ], -b[YX], -b[YY], -b[YZ], -b[ZX], -b[ZY], -b[ZZ]));
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix equality comparison
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator==(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
    const T* __RESTRICT__ b1 = m1.getBuffer();
    const T* __RESTRICT__ b2 = m2.getBuffer();
    for(int i = 0; i < 9; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix inequality operator
    @param m1 first matrix
    @param m2 second matrix */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator!=(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
    const T* __RESTRICT__ b1 = m1.getBuffer();
    const T* __RESTRICT__ b2 = m2.getBuffer();
    for(int i = 0; i < 9; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return true;
    }
    return false;
}

// -------------------------------------------------------------------------------------------------
/** @brief Matrix check for rotation
    @param m the matrix
    @param tol tolerance for numerical checks */
template <typename T>
__HOSTDEVICE__ static INLINE bool isRotation(const Matrix3<T>& m, const T tol = EPS<T>) noexcept
{
    // Check if determinant is approximately 1
    T det = determinant(m);
    if(fabs(det - T(1)) > tol)
        return false;

    // Check if matrix * transpose(matrix) = identity
    Matrix3<T> inv(transpose(m));
    inv = m * inv;

    // Check diagonal elements are approximately 1
    // clang-format off
    if(fabs(inv(XX) - T(1)) > tol || 
       fabs(inv(YY) - T(1)) > tol || 
       fabs(inv(ZZ) - T(1)) > tol)
        return false;

    // Check off-diagonal elements are approximately 0
    if(fabs(inv(XY)) > tol || fabs(inv(XZ)) > tol || fabs(inv(YX)) > tol ||
       fabs(inv(YZ)) > tol || fabs(inv(ZX)) > tol || fabs(inv(ZY)) > tol)
        return false;
    // clang-format on
    return true;
}
//@}

#endif