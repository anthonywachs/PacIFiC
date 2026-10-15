#ifndef _MISCMATH_HH_
#define _MISCMATH_HH_

#include "Basic.hh"

// =================================================================================================
/** @brief Miscellaneous math functions as header-only.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name Math functions */
//@{
/** @brief Returns whether a double number is approximately 0 with respect to a given tolerance.
    @param x the number
    @param tol the tolerance */
template <typename T>
__HOSTDEVICE__ static INLINE bool eqz(T x, T tol) noexcept
{
    return (fabs(x) <= tol);
}

// -------------------------------------------------------------------------------------------------
/** @brief Sets the minimum of 2 real numbers defined as double to these 2 numbers.
    @param x 1st real number
    @param y 2nd real number */
template <typename T>
__HOSTDEVICE__ static INLINE void set_min(T& x, T y) noexcept
{
    if(x > y)
        x = y;
}

// -------------------------------------------------------------------------------------------------
/** @brief Sets the maximum of 2 real numbers defined as double to these 2 numbers.
    @param x 1st real number
    @param y 2nd real number */
template <typename T>
__HOSTDEVICE__ static INLINE void set_max(T& x, T y) noexcept
{
    if(x < y)
        x = y;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns an angle in radians given an angle in degrees.
    @param x angle in degrees */
template <typename T>
__HOSTDEVICE__ static INLINE T rads(T x) noexcept
{
    return (x * RADS_PER_DEG<T>);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns an angle in degrees given an angle in radians.
    @param x angle in radians */
template <typename T>
__HOSTDEVICE__ static INLINE T degs(T x) noexcept
{
    return (x * DEGS_PER_RAD<T>);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the sign of the real number.
    @param x the real number */
template <typename T>
__HOSTDEVICE__ static INLINE int sgn(T x) noexcept
{
    return ((T(0) < x) - (x < T(0)));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the beta function of 2 real numbers using the gamma function, if compiling with
    NVCC, otherwise it is just a wrapper to the beta function implemented in std. The reason we
    implement this is that beta functions are not supported on device code.
    @param x the first real number
    @param y the second real number */
template <typename T>
__HOSTDEVICE__ static INLINE T grainsBeta(T x, T y) noexcept
{
#ifdef __CUDA_ARCH__
    return (tgamma(x) * tgamma(y) / tgamma(x + y));
#else
    return (std::beta(x, y));
#endif
}
//@}

#endif