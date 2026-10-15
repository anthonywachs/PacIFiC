#ifndef _QUATERNIONMATH_HH_
#define _QUATERNIONMATH_HH_

#include "Quaternion.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief Miscellaneous Quaternion functions and operators as header-only.

    Defining important quaternion functions and operators as static functions. It will increase the
    binary size, but the performance gain is much more appreciated.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name Quaternion math functions and operators */
//@{
/** @brief Returns the norm of the quaternion
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE T norm(const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    return (sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + b[3] * b[3]));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the norm squared of the quaternion
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE T norm2(const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    return (b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion conjugate
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> conjugate(const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    T __RESTRICT__        out[4];
    out[0] = -b[0];
    out[1] = -b[1];
    out[2] = -b[2];
    out[3] = b[3];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion conjugate in-place
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void conjugate(Quaternion<T>& q) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(q.getBuffer());
    b[0]              = -b[0];
    b[1]              = -b[1];
    b[2]              = -b[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion inverse
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> inverse(const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    T __RESTRICT__        out[4];
    T                     norm_inv = T(1) / norm(q);
    out[0]                         = -norm_inv * b[0];
    out[1]                         = -norm_inv * b[1];
    out[2]                         = -norm_inv * b[2];
    out[3]                         = norm_inv * b[3];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion inverse in-place
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void inverse(Quaternion<T>& q) noexcept
{
    T* __RESTRICT__ b        = const_cast<T*>(q.getBuffer());
    T               norm_inv = T(1) / norm(q);
    b[0]                     = -norm_inv * b[0];
    b[1]                     = -norm_inv * b[1];
    b[2]                     = -norm_inv * b[2];
    b[3]                     = norm_inv * b[3];
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternions addition
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator+(const Quaternion<T>& q1,
                                                     const Quaternion<T>& q2) noexcept
{
    const T* __RESTRICT__ b1 = q1.getBuffer();
    const T* __RESTRICT__ b2 = q2.getBuffer();
    T __RESTRICT__        out[4];
    for(uint i = 0; i < 4; ++i)
        out[i] = b1[i] + b2[i];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternions addition in-place
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator+=(Quaternion<T>& q1, const Quaternion<T>& q2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(q1.getBuffer());
    const T* __RESTRICT__ b2 = q2.getBuffer();
    for(uint i = 0; i < 4; ++i)
        b1[i] += b2[i];
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternions subtraction
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator-(const Quaternion<T>& q1,
                                                     const Quaternion<T>& q2) noexcept
{
    const T* __RESTRICT__ b1 = q1.getBuffer();
    const T* __RESTRICT__ b2 = q2.getBuffer();
    T __RESTRICT__        out[4];
    for(uint i = 0; i < 4; ++i)
        out[i] = b1[i] - b2[i];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternions subtraction in-place
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator-=(Quaternion<T>& q1, const Quaternion<T>& q2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(q1.getBuffer());
    const T* __RESTRICT__ b2 = q2.getBuffer();
    for(uint i = 0; i < 4; ++i)
        b1[i] -= b2[i];
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-quaternion multiplication
    @param d the multiplication factor
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator*(T d, const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    T __RESTRICT__        out[4];
    for(uint i = 0; i < 4; ++i)
        out[i] = d * b[i];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Scalar-quaternion multiplication in-place
    @param d the multiplication factor
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Quaternion<T>& q, T d) noexcept
{
    T* __RESTRICT__ b = const_cast<T*>(q.getBuffer());
    for(uint i = 0; i < 4; ++i)
        b[i] *= d;
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion-vector multiplication
    @param q the quaternion
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator*(const Quaternion<T>& q,
                                                     const Vector3<T>&    v) noexcept
{
    const T* __RESTRICT__ b1 = q.getBuffer();
    const T* __RESTRICT__ b2 = v.getBuffer();
    T __RESTRICT__        out[4];
    out[0] = (b1[3] * b2[0]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[3] * b2[1]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[3] * b2[2]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    out[3] = -(b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion-vector multiplication in-place
    @param q the quaternion
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Quaternion<T>& q, const Vector3<T>& v) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(q.getBuffer());
    const T* __RESTRICT__ b2 = v.getBuffer();
    T                     out[3];
    out[0] = (b1[3] * b2[0]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[3] * b2[1]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[3] * b2[2]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    b1[3]  = -(b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    b1[0]  = out[0];
    b1[1]  = out[1];
    b1[2]  = out[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-quaternion multiplication
    @param v the vector
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator*(const Vector3<T>&    v,
                                                     const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b1 = v.getBuffer();
    const T* __RESTRICT__ b2 = q.getBuffer();
    T __RESTRICT__        out[4];
    out[0] = (b1[0] * b2[3]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[1] * b2[3]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[2] * b2[3]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    out[3] = -(b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Vector-quaternion multiplication in-place. Note that this modifies the quaternion, not
    the vector.
    @param v the vector
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(const Vector3<T>& v, Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b1 = v.getBuffer();
    T* __RESTRICT__       b2 = const_cast<T*>(q.getBuffer());
    T                     out[3];
    out[0] = (b1[0] * b2[3]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[1] * b2[3]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[2] * b2[3]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    b2[3]  = -(b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    b2[0]  = out[0];
    b2[1]  = out[1];
    b2[2]  = out[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Rotates a vector by the inverse of a quaternion. Unit quaternion is assumed.
    @param q the quaternion
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator<<(const Quaternion<T>& q,
                                                   const Vector3<T>&    v) noexcept
{
    const T* __RESTRICT__ b1 = q.getBuffer();
    const T* __RESTRICT__ b2 = v.getBuffer();
    T __RESTRICT__        out[3];
    T                     tx = T(2) * (b1[2] * b2[1] - b1[1] * b2[2]);
    T                     ty = T(2) * (b1[0] * b2[2] - b1[2] * b2[0]);
    T                     tz = T(2) * (b1[1] * b2[0] - b1[0] * b2[1]);
    out[0]                   = b2[0] + b1[3] * tx - (b1[1] * tz - b1[2] * ty);
    out[1]                   = b2[1] + b1[3] * ty - (b1[2] * tx - b1[0] * tz);
    out[2]                   = b2[2] + b1[3] * tz - (b1[0] * ty - b1[1] * tx);
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Rotates a vector by the inverse of a quaternion in-place. Note that this modifies the
    vector, not the quaternion.
    @param v the vector
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator<<=(Vector3<T>& v, const Quaternion<T>& q) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(v.getBuffer());
    const T* __RESTRICT__ b2 = q.getBuffer();
    T                     tx = T(2) * (b2[2] * b1[1] - b2[1] * b1[2]);
    T                     ty = T(2) * (b2[0] * b1[2] - b2[2] * b1[0]);
    T                     tz = T(2) * (b2[1] * b1[0] - b2[0] * b1[1]);
    b1[0] += b2[3] * tx - (b2[1] * tz - b2[2] * ty);
    b1[1] += b2[3] * ty - (b2[2] * tx - b2[0] * tz);
    b1[2] += b2[3] * tz - (b2[0] * ty - b2[1] * tx);
}

// -------------------------------------------------------------------------------------------------
/** @brief Rotates a vector by a quaternion. Unit quaternion is assumed.
    @param q the quaternion
    @param v the vector */
template <typename T>
__HOSTDEVICE__ static INLINE Vector3<T> operator>>(const Quaternion<T>& q,
                                                   const Vector3<T>&    v) noexcept
{
    // Using the formula: v' = v + 2 * (q x v) * q + (q.w^2 - |q.v|^2) * v
    const T* __RESTRICT__ b1 = q.getBuffer();
    const T* __RESTRICT__ b2 = v.getBuffer();
    T __RESTRICT__        out[3];
    // Compute t = 2 * (q_vec x v)
    T tx = T(2) * (b1[1] * b2[2] - b1[2] * b2[1]);
    T ty = T(2) * (b1[2] * b2[0] - b1[0] * b2[2]);
    T tz = T(2) * (b1[0] * b2[1] - b1[1] * b2[0]);
    // Compute v' = v + w * t + cross(q_vec, t)
    out[0] = b2[0] + b1[3] * tx + (b1[1] * tz - b1[2] * ty);
    out[1] = b2[1] + b1[3] * ty + (b1[2] * tx - b1[0] * tz);
    out[2] = b2[2] + b1[3] * tz + (b1[0] * ty - b1[1] * tx);
    return (Vector3<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Rotates a vector by a quaternion in-place. Note that this modifies the vector, not the
    quaternion. Also, Unit quaternion is assumed.
    @param v the vector
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator>>=(Vector3<T>& v, const Quaternion<T>& q) noexcept
{
    // Using the formula: v' = v + 2 * (q x v) * q + (q.w^2 - |q.v|^2) * v
    T* __RESTRICT__       b1 = const_cast<T*>(v.getBuffer());
    const T* __RESTRICT__ b2 = q.getBuffer();
    // Compute t = 2 * (q_vec x v)
    T tx = T(2) * (b2[1] * b1[2] - b2[2] * b1[1]);
    T ty = T(2) * (b2[2] * b1[0] - b2[0] * b1[2]);
    T tz = T(2) * (b2[0] * b1[1] - b2[1] * b1[0]);
    // Compute v' = v + w * t + cross(q_vec, t)
    b1[0] += b2[3] * tx + (b2[1] * tz - b2[2] * ty);
    b1[1] += b2[3] * ty + (b2[2] * tx - b2[0] * tz);
    b1[2] += b2[3] * tz + (b2[0] * ty - b2[1] * tx);
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion-quaternion multiplication
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator*(const Quaternion<T>& q1,
                                                     const Quaternion<T>& q2) noexcept
{
    const T* __RESTRICT__ b1 = q1.getBuffer();
    const T* __RESTRICT__ b2 = q2.getBuffer();
    T __RESTRICT__        out[4];
    out[0] = (b1[3] * b2[0]) + (b1[0] * b2[3]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[3] * b2[1]) + (b1[1] * b2[3]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[3] * b2[2]) + (b1[2] * b2[3]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    out[3] = (b1[3] * b2[3]) - (b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion-quaternion multiplication in-place
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE void operator*=(Quaternion<T>& q1, const Quaternion<T>& q2) noexcept
{
    T* __RESTRICT__       b1 = const_cast<T*>(q1.getBuffer());
    const T* __RESTRICT__ b2 = q2.getBuffer();
    T                     out[3];
    out[0] = (b1[3] * b2[0]) + (b1[0] * b2[3]) + (b1[1] * b2[2]) - (b1[2] * b2[1]);
    out[1] = (b1[3] * b2[1]) + (b1[1] * b2[3]) + (b1[2] * b2[0]) - (b1[0] * b2[2]);
    out[2] = (b1[3] * b2[2]) + (b1[2] * b2[3]) + (b1[0] * b2[1]) - (b1[1] * b2[0]);
    b1[3]  = (b1[3] * b2[3]) - (b1[0] * b2[0]) - (b1[1] * b2[1]) - (b1[2] * b2[2]);
    b1[0]  = out[0];
    b1[1]  = out[1];
    b1[2]  = out[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion equality operator
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator==(const Quaternion<T>& q1,
                                             const Quaternion<T>& q2) noexcept
{
    const T* __RESTRICT__ b1 = q1.getBuffer();
    const T* __RESTRICT__ b2 = q2.getBuffer();
    for(int i = 0; i < 4; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion inequality operator
    @param q1 1st quaternion
    @param q2 2nd quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE bool operator!=(const Quaternion<T>& q1,
                                             const Quaternion<T>& q2) noexcept
{
    const T* __RESTRICT__ b1 = q1.getBuffer();
    const T* __RESTRICT__ b2 = q2.getBuffer();
    for(int i = 0; i < 4; ++i)
    {
        if(fabs(b1[i] - b2[i]) > EPS<T>)
            return true;
    }
    return false;
}

// -------------------------------------------------------------------------------------------------
/** @brief Quaternion sign flip
    @param q the quaternion */
template <typename T>
__HOSTDEVICE__ static INLINE Quaternion<T> operator-(const Quaternion<T>& q) noexcept
{
    const T* __RESTRICT__ b = q.getBuffer();
    T __RESTRICT__        out[4];
    for(uint i = 0; i < 4; ++i)
        out[i] = -b[i];
    return (Quaternion<T>(out));
}

// -------------------------------------------------------------------------------------------------
/** @brief Fused Math for Minkowski difference of two points; w = a - q_b2a(b) + v.
    @param a the first point
    @param b the second point
    @param v the relative position vector from body2 to body1
    @param q the quaternion representing the rotation from body2 to body1 */
template <typename T>
__HOSTDEVICE__ static INLINE void FusedMinkowskiDifference(const Vector3<T>&    a,
                                                           const Vector3<T>&    b,
                                                           const Vector3<T>&    v,
                                                           const Quaternion<T>& q,
                                                           Vector3<T>&          w) noexcept
{
    const T* __RESTRICT__ bb = b.getBuffer();
    const T* __RESTRICT__ bq = q.getBuffer();
    T                     tx = T(2) * (bq[1] * bb[2] - bq[2] * bb[1]);
    T                     ty = T(2) * (bq[2] * bb[0] - bq[0] * bb[2]);
    T                     tz = T(2) * (bq[0] * bb[1] - bq[1] * bb[0]);

    const T* __RESTRICT__ ba = a.getBuffer();
    const T* __RESTRICT__ bv = v.getBuffer();
    T* __RESTRICT__       bw = const_cast<T*>(w.getBuffer());
    bw[0]                    = ba[0] - bb[0] - bq[3] * tx - bq[1] * tz + bq[2] * ty - bv[0];
    bw[1]                    = ba[1] - bb[1] - bq[3] * ty - bq[2] * tx + bq[0] * tz - bv[1];
    bw[2]                    = ba[2] - bb[2] - bq[3] * tz - bq[0] * ty + bq[1] * tx - bv[2];
}

// -------------------------------------------------------------------------------------------------
/** @brief Fused Math for Minkowski difference of two points;
    w = q_a2w(a) - q_b2w(b) + v_a2w - v_b2w.
    @param a the first point
    @param b the second point
    @param v_a2w the position vector of body1 in world frame
    @param v_b2w the position vector of body2 in world frame
    @param q_a2w the quaternion representing body1 to world
    @param q_b2w the quaternion representing body2 to world
    @param w the resulting Minkowski difference */
template <typename T>
__HOSTDEVICE__ static INLINE void FusedMinkowskiDifference(const Vector3<T>&    a,
                                                           const Vector3<T>&    b,
                                                           const Vector3<T>&    v_a2w,
                                                           const Vector3<T>&    v_b2w,
                                                           const Quaternion<T>& q_a2w,
                                                           const Quaternion<T>& q_b2w,
                                                           Vector3<T>&          w) noexcept
{
    // w = (q_a2w >> a) - (q_b2w >> b) + v_a2w - v_b2w;
    const T* __RESTRICT__ ba  = a.getBuffer();
    const T* __RESTRICT__ bqa = q_a2w.getBuffer();
    T                     txa = T(2) * (bqa[1] * ba[2] - bqa[2] * ba[1]);
    T                     tya = T(2) * (bqa[2] * ba[0] - bqa[0] * ba[2]);
    T                     tza = T(2) * (bqa[0] * ba[1] - bqa[1] * ba[0]);

    const T* __RESTRICT__ bb  = b.getBuffer();
    const T* __RESTRICT__ bqb = q_b2w.getBuffer();
    T                     txb = T(2) * (bqb[1] * bb[2] - bqb[2] * bb[1]);
    T                     tyb = T(2) * (bqb[2] * bb[0] - bqb[0] * bb[2]);
    T                     tzb = T(2) * (bqb[0] * bb[1] - bqb[1] * bb[0]);

    const T* __RESTRICT__ bva = v_a2w.getBuffer();
    const T* __RESTRICT__ bvb = v_b2w.getBuffer();
    T* __RESTRICT__       bw  = const_cast<T*>(w.getBuffer());
    // clang format off
    bw[0] = ba[0] + bqa[3] * txa + bqa[1] * tza - bqa[2] * tya + bva[0] - bb[0] - bqb[3] * txb
            - bqb[1] * tzb + bqb[2] * tyb - bvb[0];
    bw[1] = ba[1] + bqa[3] * tya + bqa[2] * txa - bqa[0] * tza + bva[1] - bb[1] - bqb[3] * tyb
            - bqb[2] * txb + bqb[0] * tzb - bvb[1];
    bw[2] = ba[2] + bqa[3] * tza + bqa[0] * tya - bqa[1] * txa + bva[2] - bb[2] - bqb[3] * tzb
            - bqb[0] * tyb + bqb[1] * txb - bvb[2];
    // clang format on
}

// -------------------------------------------------------------------------------------------------
/** @brief Transforms a point using a quaternion and a vector.
    @param q the quaternion representing the rotation
    @param v the translation vector
    @param w the vector to be transformed */
template <typename T>
__HOSTDEVICE__ static INLINE void
    transform(const Quaternion<T>& q, const Vector3<T>& v, Vector3<T>& w) noexcept
{
    const T* __RESTRICT__ bq = q.getBuffer();
    T* __RESTRICT__       bw = const_cast<T*>(w.getBuffer());
    T                     tx = T(2) * (bq[1] * bw[2] - bq[2] * bw[1]);
    T                     ty = T(2) * (bq[2] * bw[0] - bq[0] * bw[2]);
    T                     tz = T(2) * (bq[0] * bw[1] - bq[1] * bw[0]);

    const T* __RESTRICT__ bv = v.getBuffer();
    bw[0] += bq[3] * tx + (bq[1] * tz - bq[2] * ty) + bv[0];
    bw[1] += bq[3] * ty + (bq[2] * tx - bq[0] * tz) + bv[1];
    bw[2] += bq[3] * tz + (bq[0] * ty - bq[1] * tx) + bv[2];
}
//@}

#endif