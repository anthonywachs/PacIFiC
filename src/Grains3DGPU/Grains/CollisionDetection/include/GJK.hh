#ifndef _GJK_JH_HH_
#define _GJK_JH_HH_

#include "Convex.hh"
#include "GJK_ShapeData.hh"
#include "Quaternion.hh"
#include "Transform3.hh"

/** @brief GJK sub-algorithm variant */
enum class GJKType
{
    JOHNSON,     /**< Johnson's algorithm */
    SIGNEDVOLUME /**< Signed volume algorithm */
};

// =================================================================================================
/** @brief The header for the GJK distance query algorithm.

    The GJK distance query algorithm using the Johnson and signedVolume subalgorithm with the
    backup procedure. It supports both single and double floating point operations, but it is not
    recommended to use the single precision version as it is prone to numerical instabilities.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name GJK : External methods */
//@{
/** @brief Returns whether 2 convex shapes intersect - relative transformation.
    @param a convex shape A
    @param b convex shape B
    @param b2a geometric transformation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>& a, const Convex<T>& b, const Transform3<T>& b2a);

/** @brief Returns whether 2 convex shapes intersect.
    @param a convex shape A
    @param b convex shape B
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Transform3<T>& a2w,
                                 const Transform3<T>& b2w);

/** @brief Returns whether 2 convex shapes intersect - relative transformation.
    @param a convex shape A
    @param b convex shape B
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Vector3<T>&    v_b2a,
                                 const Quaternion<T>& q_b2a);

/** @brief Returns whether 2 convex shapes intersect.
    @param a convex shape A
    @param b convex shape B
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Vector3<T>&    v_a2w,
                                 const Vector3<T>&    v_b2w,
                                 const Quaternion<T>& q_a2w,
                                 const Quaternion<T>& q_b2w);

/** @brief Returns the minimal distance between 2 convex shapes and a point per convex shape that
    represents the tips of the minimal distance segment -- relative transformation.
    @param a convex shape A
    @param b convex shape B
    @param b2a geometric transformation describing convex B in the reference frame of A
    @param crustA crust/skin thickness on A (shrinks A along search dir)
    @param crustB crust/skin thickness on B (shrinks B along search dir)
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the other tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Transform3<T>& b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);

/** @brief Returns the minimal distance between 2 convex shapes and a point per convex shape that
    represents the tips of the minimal distance segment.
    @param a convex shape A
    @param b convex shape B
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame
    @param crustA crust/skin thickness on A (shrinks A along search dir)
    @param crustB crust/skin thickness on B (shrinks B along search dir)
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Transform3<T>& a2w,
                                          const Transform3<T>& b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);

/** @brief Returns the minimal distance between 2 convex shapes and a point per convex shape that
    represents the tips of the minimal distance segment -- relative transformation.
    @param a convex shape A
    @param b convex shape B
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param crustA crust/skin thickness on A (shrinks A along search dir)
    @param crustB crust/skin thickness on B (shrinks B along search dir)
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Vector3<T>&    v_b2a,
                                          const Quaternion<T>& q_b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);

/** @brief Returns the minimal distance between 2 convex shapes and a point per convex shape that
    represents the tips of the minimal distance segment.
    @param a convex shape A
    @param b convex shape B
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param crustA crust/skin thickness on A (shrinks A along search dir)
    @param crustB crust/skin thickness on B (shrinks B along search dir)
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Vector3<T>&    v_a2w,
                                          const Vector3<T>&    v_b2w,
                                          const Quaternion<T>& q_a2w,
                                          const Quaternion<T>& q_b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);

/** @brief Returns the minimal distance between 2 convex shapes using pre-built ShapeData for
    vtable-free support evaluation -- relative transformation (vec/quat).
    @param sdA ShapeData for shape A
    @param sdB ShapeData for shape B
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param crustA crust thickness on A
    @param crustB crust thickness on B
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const ShapeData<T>&  sdA,
                                          const ShapeData<T>&  sdB,
                                          const Vector3<T>&    v_b2a,
                                          const Quaternion<T>& q_b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);

/** @brief Returns the minimal distance between 2 convex shapes using pre-built ShapeData for
    vtable-free support evaluation -- world-frame (vec/quat).
    @param sdA ShapeData for shape A
    @param sdB ShapeData for shape B
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param crustA crust thickness on A
    @param crustB crust thickness on B
    @param pa point representing one tip of the minimal distance segment on A
    @param pb point representing the tip of the minimal distance segment on B
    @param nbIter number of iterations of GJK for convergence */
template <typename T, GJKType GJKType, bool Acceleration = false>
__HOSTDEVICE__ T computeClosestPoints_GJK(const ShapeData<T>&  sdA,
                                          const ShapeData<T>&  sdB,
                                          const Vector3<T>&    v_a2w,
                                          const Vector3<T>&    v_b2w,
                                          const Quaternion<T>& q_a2w,
                                          const Quaternion<T>& q_b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter);
//@}

#endif