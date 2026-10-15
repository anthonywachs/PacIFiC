#ifndef _OBB_HH_
#define _OBB_HH_

#include "BoundingBox.hh"
#include "MatrixMath.hh"
#include "Quaternion.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The header for the axis-aligned and oriented bounding boxes collision detection.

    Axis-aligned Bounding Boxes (AABB) and Oriented Bounding Boxes (OBB) routines to find whether
    bounding boxes are in contact or not. AABB is deprecated, so try to use OBB.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name OBB: External methods */
//@{
// Low-level methods for OBB as macros in double precision
#define TESTCASE1(i) \
    (fabs(cen[i]) > (a[i] + b[0] * oriAbs[i][0] + b[1] * oriAbs[i][1] + b[2] * oriAbs[i][2]))

#define TESTCASE2(i)                                                    \
    (fabs(cen[0] * ori[0][i] + cen[1] * ori[1][i] + cen[2] * ori[2][i]) \
     > (b[i] + a[0] * oriAbs[0][i] + a[1] * oriAbs[1][i] + a[2] * oriAbs[2][i]))

#define TESTCASE3(i, j)                                                                    \
    (fabs(cen[(i + 2) % 3] * ori[(i + 1) % 3][j] - cen[(i + 1) % 3] * ori[(i + 2) % 3][j]) \
     > (a[(i + 1) % 3] * oriAbs[(i + 2) % 3][j] + a[(i + 2) % 3] * oriAbs[(i + 1) % 3][j]  \
        + b[(j + 1) % 3] * oriAbs[i][(j + 2) % 3] + b[(j + 2) % 3] * oriAbs[i][(j + 1) % 3]))

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using OBB test - absolute
    transformation.
    @param a The first bounding box
    @param b The second bounding box
    @param trA2W The transformation from A's local space to world space
    @param trB2W The transformation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingBox(const Vector3<T>&    a,
                                                 const Vector3<T>&    b,
                                                 const Transform3<T>& trA2W,
                                                 const Transform3<T>& trB2W)
{
    // First, we compute the transpose of trA2W basis and store it in ori
    Matrix3<T> ori(transpose(trA2W.getBasis()));
    // Then, the center is
    const Vector3<T>& cen = ori * (trB2W.getOrigin() - trA2W.getOrigin());
    // Finally, we compute the actual relative rotation matrix
    ori *= trB2W.getBasis();
    // And, we compute the absolute value of the matrix + some noise to
    // encounter arithmetic errors.
    const Matrix3<T> oriAbs(fabs(ori(0, 0)) + LOWEPS<T>,
                            fabs(ori(0, 1)) + LOWEPS<T>,
                            fabs(ori(0, 2)) + LOWEPS<T>,
                            fabs(ori(1, 0)) + LOWEPS<T>,
                            fabs(ori(1, 1)) + LOWEPS<T>,
                            fabs(ori(1, 2)) + LOWEPS<T>,
                            fabs(ori(2, 0)) + LOWEPS<T>,
                            fabs(ori(2, 1)) + LOWEPS<T>,
                            fabs(ori(2, 2)) + LOWEPS<T>);

    // CASE 1: ( three of them )
    if TESTCASE1(0)
        return (false);
    if TESTCASE1(1)
        return (false);
    if TESTCASE1(2)
        return (false);

    // CASE 2: ( three of them )
    if TESTCASE2(0)
        return (false);
    if TESTCASE2(1)
        return (false);
    if TESTCASE2(2)
        return (false);

    // CASE 3: ( nine of them )
    if TESTCASE3(0, 0)
        return (false);
    if TESTCASE3(1, 0)
        return (false);
    if TESTCASE3(2, 0)
        return (false);
    if TESTCASE3(0, 1)
        return (false);
    if TESTCASE3(1, 1)
        return (false);
    if TESTCASE3(2, 1)
        return (false);
    if TESTCASE3(0, 2)
        return (false);
    if TESTCASE3(1, 2)
        return (false);
    if TESTCASE3(2, 2)
        return (false);

    return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using OBB test - relative
    transformation.
    @param a The first bounding box
    @param b The second bounding box
    @param trB2A The transformation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingBox(const Vector3<T>&    a,
                                                 const Vector3<T>&    b,
                                                 const Transform3<T>& trB2A)
{
    const Vector3<T>& cen = trB2A.getOrigin();
    const Matrix3<T>& ori = trB2A.getBasis();
    Matrix3<T> const  oriAbs(fabs(ori(0, 0)) + LOWEPS<T>,
                            fabs(ori(0, 1)) + LOWEPS<T>,
                            fabs(ori(0, 2)) + LOWEPS<T>,
                            fabs(ori(1, 0)) + LOWEPS<T>,
                            fabs(ori(1, 1)) + LOWEPS<T>,
                            fabs(ori(1, 2)) + LOWEPS<T>,
                            fabs(ori(2, 0)) + LOWEPS<T>,
                            fabs(ori(2, 1)) + LOWEPS<T>,
                            fabs(ori(2, 2)) + LOWEPS<T>);

    // CASE 1: ( three of them )
    if TESTCASE1(0)
        return (false);
    if TESTCASE1(1)
        return (false);
    if TESTCASE1(2)
        return (false);

    // CASE 2: ( three of them )
    if TESTCASE2(0)
        return (false);
    if TESTCASE2(1)
        return (false);
    if TESTCASE2(2)
        return (false);

    // CASE 3: ( nine of them )
    if TESTCASE3(0, 0)
        return (false);
    if TESTCASE3(1, 0)
        return (false);
    if TESTCASE3(2, 0)
        return (false);
    if TESTCASE3(0, 1)
        return (false);
    if TESTCASE3(1, 1)
        return (false);
    if TESTCASE3(2, 1)
        return (false);
    if TESTCASE3(0, 2)
        return (false);
    if TESTCASE3(1, 2)
        return (false);
    if TESTCASE3(2, 2)
        return (false);

    return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using OBB test - quaternion version.
    @param a The first bounding box
    @param b The second bounding box
    @param v_a2w The translation from A's local space to world space
    @param v_b2w The translation from B's local space to world space
    @param q_a2w The rotation from A's local space to world space
    @param q_b2w The rotation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingBox(const Vector3<T>&    a,
                                                 const Vector3<T>&    b,
                                                 const Vector3<T>&    v_a2w,
                                                 const Vector3<T>&    v_b2w,
                                                 const Quaternion<T>& q_a2w,
                                                 const Quaternion<T>& q_b2w)
{
    // Compute relative rotation: q_b2a = inverse(q_a2w) * q_b2w
    Matrix3<T> ori = (inverse(q_a2w) * q_b2w).toMatrix();
    // Relative translation: v_b2a = q_a2w << (v_b2w - v_a2w)
    Vector3<T> cen = q_a2w << (v_b2w - v_a2w);

    // Compute absolute value matrix
    const Matrix3<T> oriAbs(fabs(ori(0, 0)) + LOWEPS<T>,
                            fabs(ori(0, 1)) + LOWEPS<T>,
                            fabs(ori(0, 2)) + LOWEPS<T>,
                            fabs(ori(1, 0)) + LOWEPS<T>,
                            fabs(ori(1, 1)) + LOWEPS<T>,
                            fabs(ori(1, 2)) + LOWEPS<T>,
                            fabs(ori(2, 0)) + LOWEPS<T>,
                            fabs(ori(2, 1)) + LOWEPS<T>,
                            fabs(ori(2, 2)) + LOWEPS<T>);

    // CASE 1: ( three of them )
    if TESTCASE1(0)
        return (false);
    if TESTCASE1(1)
        return (false);
    if TESTCASE1(2)
        return (false);

    // CASE 2: ( three of them )
    if TESTCASE2(0)
        return (false);
    if TESTCASE2(1)
        return (false);
    if TESTCASE2(2)
        return (false);

    // CASE 3: ( nine of them )
    if TESTCASE3(0, 0)
        return (false);
    if TESTCASE3(1, 0)
        return (false);
    if TESTCASE3(2, 0)
        return (false);
    if TESTCASE3(0, 1)
        return (false);
    if TESTCASE3(1, 1)
        return (false);
    if TESTCASE3(2, 1)
        return (false);
    if TESTCASE3(0, 2)
        return (false);
    if TESTCASE3(1, 2)
        return (false);
    if TESTCASE3(2, 2)
        return (false);

    return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using OBB test - quaternion relative
    version.
    @param a The first bounding box
    @param b The second bounding box
    @param v_b2a The translation from B's local space to A's local space
    @param q_b2a The rotation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingBox(const Vector3<T>&    a,
                                                 const Vector3<T>&    b,
                                                 const Vector3<T>&    v_b2a,
                                                 const Quaternion<T>& q_b2a)
{
    const Vector3<T>& cen = v_b2a;
    Matrix3<T>        ori = q_b2a.toMatrix();

    Matrix3<T> const oriAbs(fabs(ori(0, 0)) + LOWEPS<T>,
                            fabs(ori(0, 1)) + LOWEPS<T>,
                            fabs(ori(0, 2)) + LOWEPS<T>,
                            fabs(ori(1, 0)) + LOWEPS<T>,
                            fabs(ori(1, 1)) + LOWEPS<T>,
                            fabs(ori(1, 2)) + LOWEPS<T>,
                            fabs(ori(2, 0)) + LOWEPS<T>,
                            fabs(ori(2, 1)) + LOWEPS<T>,
                            fabs(ori(2, 2)) + LOWEPS<T>);

    // CASE 1: ( three of them )
    if TESTCASE1(0)
        return (false);
    if TESTCASE1(1)
        return (false);
    if TESTCASE1(2)
        return (false);

    // CASE 2: ( three of them )
    if TESTCASE2(0)
        return (false);
    if TESTCASE2(1)
        return (false);
    if TESTCASE2(2)
        return (false);

    // CASE 3: ( nine of them )
    if TESTCASE3(0, 0)
        return (false);
    if TESTCASE3(1, 0)
        return (false);
    if TESTCASE3(2, 0)
        return (false);
    if TESTCASE3(0, 1)
        return (false);
    if TESTCASE3(1, 1)
        return (false);
    if TESTCASE3(2, 1)
        return (false);
    if TESTCASE3(0, 2)
        return (false);
    if TESTCASE3(1, 2)
        return (false);
    if TESTCASE3(2, 2)
        return (false);

    return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using AABB test - absolute
    transformation.
    @param a The first bounding box
    @param b The second bounding box
    @param trA2W The transformation from A's local space to world space
    @param trB2W The transformation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectAxisAlignedBoundingBox(const Vector3<T>&    a,
                                                    const Vector3<T>&    b,
                                                    const Transform3<T>& trA2W,
                                                    const Transform3<T>& trB2W)
{
    // TODO: a and b should be modified according to trA2W and trB2W
    // TODO: should we do len = a.getExtent() + b.getExtent()?
    const Vector3<T>& posA = trA2W.getOrigin();
    const Vector3<T>& posB = trB2W.getOrigin();
    if(fabs(posA[X] - posB[X]) > (a[X] + b[X]))
        return (false);
    else if(fabs(posA[Y] - posB[Y]) > (a[Y] + b[Y]))
        return (false);
    else if(fabs(posA[Z] - posB[Z]) > (a[Z] + b[Z]))
        return (false);
    else  // overlap
        return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using AABB test - relative
    transformation.
    @param a The first bounding box
    @param b The second bounding box
    @param trB2A The transformation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectAxisAlignedBoundingBox(const Vector3<T>&    a,
                                                    const Vector3<T>&    b,
                                                    const Transform3<T>& trB2A)
{
    // TODO: a and b should be modified according to trA2W and trB2W
    // TODO: should we do len = a.getExtent() + b.getExtent()?
    const Vector3<T>& pos = trB2A.getOrigin();
    if(fabs(pos[X]) > (a[X] + b[X]))
        return (false);
    else if(fabs(pos[Y]) > (a[Y] + b[Y]))
        return (false);
    else if(fabs(pos[Z]) > (a[Z] + b[Z]))
        return (false);
    else  // overlap
        return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using AABB test - quaternion version.
    @param a The first bounding box
    @param b The second bounding box
    @param v_a2w The translation from A's local space to world space
    @param v_b2w The translation from B's local space to world space
    @param q_a2w The rotation from A's local space to world space
    @param q_b2w The rotation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectAxisAlignedBoundingBox(const Vector3<T>&    a,
                                                    const Vector3<T>&    b,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w)
{
    // For AABB, we ignore rotation and just check axis-aligned extents
    if(fabs(v_a2w[X] - v_b2w[X]) > (a[X] + b[X]))
        return (false);
    else if(fabs(v_a2w[Y] - v_b2w[Y]) > (a[Y] + b[Y]))
        return (false);
    else if(fabs(v_a2w[Z] - v_b2w[Z]) > (a[Z] + b[Z]))
        return (false);
    else  // overlap
        return (true);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether the bounding boxes are in contact using AABB test - quaternion relative
    version.
    @param a The first bounding box
    @param b The second bounding box
    @param v_b2a The translation from B's local space to A's local space
    @param q_b2a The rotation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectAxisAlignedBoundingBox(const Vector3<T>&    a,
                                                    const Vector3<T>&    b,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a)
{
    // For AABB, we ignore rotation and just check axis-aligned extents
    if(fabs(v_b2a[X]) > (a[X] + b[X]))
        return (false);
    else if(fabs(v_b2a[Y]) > (a[Y] + b[Y]))
        return (false);
    else if(fabs(v_b2a[Z]) > (a[Z] + b[Z]))
        return (false);
    else  // overlap
        return (true);
}

// Undefining the low-level methods
#undef TESTCASE1
#undef TESTCASE2
#undef TESTCASE3
//@}

#endif