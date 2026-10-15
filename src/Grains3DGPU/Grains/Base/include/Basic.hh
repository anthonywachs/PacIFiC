#ifndef _BASIC_HH_
#define _BASIC_HH_

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <algorithm>
#include <bits/stdc++.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <execution>
#include <float.h>
#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <new>
#include <ostream>
#include <random>
#include <stdio.h>
#include <string>

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include <omp.h>

#include "ReaderXML.hh"
#include "WriterXML.hh"

// =================================================================================================
/** @brief Various constants and type definitions.

    @author A.Yazdani - 2023 - Construction */
// =================================================================================================
/** @name Macros */
//@{
/** @brief Compiler macros */
#if defined(__NVCC__) || defined(__CUDACC__)
#define __HOST__ __host__
#define __DEVICE__ __device__
#define __HOSTDEVICE__ __host__ __device__
#define __MANAGED__ __managed__
#define __GLOBAL__ __global__
#define INLINE __inline__
#define __RESTRICT__ __restrict__
#else
// For testing or non-CUDA compilation, use simplified macros
#define __HOST__
#define __DEVICE__
#define __HOSTDEVICE__
#define __GLOBAL__
#define INLINE inline
#define __RESTRICT__
#endif

// -------------------------------------------------------------------------------------------------
/** @brief Type macros */
//@}

/** @name Typedefs */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Unsigned integer type */
using uint = unsigned int;
//@}

/** @name Enumerations */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Space dimensions */
enum Direction
{
    X,    // x direction
    Y,    // y direction
    Z,    // z direction
    W,    // scalar component of quaternions
    NONE  // no direction
};

// -------------------------------------------------------------------------------------------------
/** @brief Matrix dimensions */
enum MatDirection
{
    XX,
    XY,
    XZ,
    YX,
    YY,
    YZ,
    ZX,
    ZY,
    ZZ
};
//@}

/** @name Constant expressions at compile-time */
//@{
// -------------------------------------------------------------------------------------------------
/** \brief PI value */
template <class T>
constexpr T PI = T(3.1415926535897932385L);

// -------------------------------------------------------------------------------------------------
/** \brief 2*PI value */
template <class T>
constexpr T TWO_PI = T(6.28318530717958623200L);

// -------------------------------------------------------------------------------------------------
/** \brief Degree per radian value */
template <class T>
constexpr T DEGS_PER_RAD = T(57.29577951308232286465L);

// -------------------------------------------------------------------------------------------------
/** \brief Radian per degree value */
template <class T>
constexpr T RADS_PER_DEG = T(0.01745329251994329547L);

// -------------------------------------------------------------------------------------------------
/** \brief High (Machine) epsilon value */
template <class T>
constexpr T HIGHEPS = T(1.e-15);
template <>
constexpr float HIGHEPS<float> = float(1.e-08);

// -------------------------------------------------------------------------------------------------
/** \brief Epsilon value */
template <class T>
constexpr T EPS = T(1.e-10);
template <>
constexpr float EPS<float> = float(1.e-05);

// -------------------------------------------------------------------------------------------------
/** \brief Low epsilon value */
template <class T>
constexpr T LOWEPS = T(1.e-05);
template <>
constexpr float LOWEPS<float> = float(1.e-03);
//@}

#endif