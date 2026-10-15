#ifndef _COLLISIONDETECTION_HH_
#define _COLLISIONDETECTION_HH_

#include <tuple>

#include "ContactInfo.hh"
#include "GJK.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "MatrixMath.hh"
#include "MiscMath.hh"
#include "OBB.hh"
#include "QuaternionMath.hh"
#include "Rectangle.hh"
#include "RigidBody.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The header-only file for Rigid bodies collision detections.

    Functions for collision detection between two rigid bodies.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name CollisionDetection Low-Level Methods */
//@{
/** @brief Helper that extracts commonly-used rigid-body properties for contact handling.
    @param rbA first rigid body
    @param rbB second rigid body */
template <typename T>
__HOSTDEVICE__ static INLINE auto getPropertiesForContact(const RigidBody<T>& rbA,
                                                          const RigidBody<T>& rbB)
{
    auto sa = rbA.getPropertiesSnapshot();
    auto sb = rbB.getPropertiesSnapshot();

    const Convex<T>* convexA     = sa.convex;
    const Convex<T>* convexB     = sb.convex;
    T                crustA      = sa.crustThickness;
    T                crustB      = sb.crustThickness;
    T                circRadiusA = sa.circumscribedRadius;
    T                circRadiusB = sb.circumscribedRadius;

    T massA          = sa.mass;
    T massB          = sb.mass;
    T invMassA       = (massA == T(0)) ? T(0) : T(1) / massA;
    T invMassB       = (massB == T(0)) ? T(0) : T(1) / massB;
    T invReducedMass = 1 / (invMassA + invMassB);

    T    invRadA       = (circRadiusA == T(0)) ? T(0) : T(1) / circRadiusA;
    T    invRadB       = (circRadiusB == T(0)) ? T(0) : T(1) / circRadiusB;
    T    averageRadius = 1 / (invRadA + invRadB);
    uint materialHash  = triangularHash(sa.material, sb.material);

    return std::make_tuple(convexA,
                           crustA,
                           circRadiusA,
                           convexB,
                           crustB,
                           circRadiusB,
                           invReducedMass,
                           averageRadius,
                           materialHash);
}

// -------------------------------------------------------------------------------------------------
/** @brief ShapeData-only variant of getPropertiesForContact. Reads all needed scalar data
    directly from ShapeData fields -- no RigidBody or virtual dispatch required.
    @param sdA pre-built ShapeData for body A
    @param sdB pre-built ShapeData for body B */
template <typename T>
__HOSTDEVICE__ static INLINE auto getPropertiesForContact(const ShapeData<T>& sdA,
                                                          const ShapeData<T>& sdB)
{
    const T    crustA         = sdA.crust;
    const T    crustB         = sdB.crust;
    const T    circRadiusA    = sdA.circumscribedRadius;
    const T    circRadiusB    = sdB.circumscribedRadius;
    const T    invReducedMass = T(1) / (sdA.invMass + sdB.invMass);
    const T    invRadA        = (circRadiusA == T(0)) ? T(0) : T(1) / circRadiusA;
    const T    invRadB        = (circRadiusB == T(0)) ? T(0) : T(1) / circRadiusB;
    const T    averageRadius  = T(1) / (invRadA + invRadB);
    const uint materialHash   = triangularHash(sdA.material, sdB.material);
    return std::make_tuple(crustA,
                           circRadiusA,
                           crustB,
                           circRadiusB,
                           invReducedMass,
                           averageRadius,
                           materialHash);
}

// -------------------------------------------------------------------------------------------------
// Returns whether two rigid bodies os spherical shape intersect
template <typename T>
__HOSTDEVICE__ static INLINE bool
    intersectSpheres(const RigidBody<T>& rbA, const RigidBody<T>& rbB, const Vector3<T>& v_b2a)
{
    T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
    T dist2    = norm2(v_b2a);
    return (dist2 < radiiSum * radiiSum);
}
//@}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
/** @name CollisionDetection High-Level Methods */
//@{
/** @brief Returns whether 2 rigid bodies intersect - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param b2a geometric transformation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ inline bool
    intersectRigidBodies(const RigidBody<T>& rbA, const RigidBody<T>& rbB, const Transform3<T>& b2a)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, b2a));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Transform3<T>& a2w,
                                                const Transform3<T>& b2w)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, a2w, b2w));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Vector3<T>&    v_b2a,
                                                const Quaternion<T>& q_b2a)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, v_b2a, q_b2a));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether 2 rigid bodies intersect.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame */
template <typename T>
__HOSTDEVICE__ inline bool intersectRigidBodies(const RigidBody<T>&  rbA,
                                                const RigidBody<T>&  rbB,
                                                const Vector3<T>&    v_a2w,
                                                const Vector3<T>&    v_b2w,
                                                const Quaternion<T>& q_a2w,
                                                const Quaternion<T>& q_b2w)
{
    const Convex<T>& convexA = *(rbA.getConvex());
    const Convex<T>& convexB = *(rbB.getConvex());
    return (intersectGJK(convexA, convexB, v_a2w, v_b2w, q_a2w, q_b2w));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param b2a geometric transformation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Transform3<T>& b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2a.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    /* ---------------------------------------------------------------------------------------------
    Comments on the contactInfo. It applies to all variants of this function:
    1. If actual overlap distance (GJK dist - crustA - crustB < 0), there is contact otherwise no
    contact. Although we can enforce an early exit, but other threads are most likely still
    running, so we continue to have consistent code path.

    2. ptA and ptB are in their respective local coordinate systems and represent points on the
    actual rigid bodies, not the shrunken versions. Contact point definition as the mid point
    between ptA and ptB.

    3. If contact, overlap is negative and overlap_vector is from B to A If no contact, overlap is
    positive and we do not care about the direction of overlap_vector. Assuming A and B are the
    centers of the 2 convex bodies overlap_vector = overlap * Vector3(A to B)
    --------------------------------------------------------------------------------------------- */
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = b2a.getOrigin();
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Rectangle-Particle Case
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        Vector3<T>       ptA, ptB;
        const Vector3<T> r = b2a.getOrigin()[Z] > 0 ? Vector3<T>(0, 0, -1) : Vector3<T>(0, 0, 1);
        ptB                = (b2a)(convexB.support(r * b2a.getBasis()));
        if(ptB[Z] < T(0))
        {
            ptA = Vector3<T>(ptB[X], ptB[Y], T(0));
            if(convexA.isInside(ptA))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case
    else
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        // ptA = (a2a)(ptA);
        ptB                    = (b2a)(ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Transform3<T>& a2w,
                                                    const Transform3<T>& b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2w.getOrigin() - a2w.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         a2w,
                                         b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(b2w.getOrigin() - a2w.getOrigin()) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              a2w,
                                              b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = a2w.getOrigin();
        ptB                      = b2w.getOrigin() - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Rectangle-Particle Case
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        const Vector3<T>& c = a2w.getOrigin();
        const Matrix3<T>& m = a2w.getBasis();
        // rectangle normal is a2w.getBasis() * [0, 0, 1] which is the last column of the transform
        Vector3<T> r(m(XZ), m(YZ), m(ZZ));
        r.normalize();
        r *= copysign(T(1), r * (b2w.getOrigin() - c));
        ptB = (b2w)(convexB.support((-r) * b2w.getBasis()));
        if(r * (ptB - c) < T(0))
        {
            ptA = ((c - ptB) * r) * r + ptB;
            if(convexA.isInside(inverse(m) * (ptA - c)))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case
    else
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   a2w,
                                                                                   b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        ptA                    = (a2w)(ptA);
        ptB                    = (b2w)(ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies - relative transformation.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2a) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_b2a,
                                         q_b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Rectangle-Particle Case
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        Vector3<T>       ptA, ptB;
        const Vector3<T> r = v_b2a[Z] > 0 ? Vector3<T>(0, 0, -1) : Vector3<T>(0, 0, 1);
        ptB                = convexB.support(q_b2a << r);
        transform(q_b2a, v_b2a, ptB);
        if(ptB[Z] < T(0))
        {
            ptA = Vector3<T>(ptB[X], ptB[Y], T(0));
            if(convexA.isInside(ptA))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case
    else
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        // transform(q_a2a, v_a2a, ptA);
        transform(q_b2a, v_b2a, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_a2w,
                                         v_b2w,
                                         q_a2w,
                                         q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_a2w,
                                              v_b2w,
                                              q_a2w,
                                              q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Rectangle-Particle Case
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        Vector3<T> r = q_a2w >> Vector3<T>(0, 0, 1);
        r.normalize();
        r *= copysign(T(1), r * (v_b2w - v_a2w));
        ptB = (q_b2w >> convexB.support(q_b2w << (-r))) + v_b2w;
        if(r * (ptB - v_a2w) < T(0))
        {
            ptA = ((v_a2w - ptB) * r) * r + ptB;
            if(convexA.isInside(q_a2w << (ptA - v_a2w)))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case
    else
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(convexA,
                                                                                   convexB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        transform(q_a2w, v_a2w, ptA);
        transform(q_b2w, v_b2w, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the distance between 2 rigid bodies.
    @param rbA first rigid body
    @param rbB second rigid body
    @param a2w geometric transformation describing convex A in the world reference frame
    @param b2w geometric transformation describing convex B in the world reference frame
    @param method method identifier (currently unused) */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline T distanceRigidBodies(const RigidBody<T>&  rbA,
                                            const RigidBody<T>&  rbB,
                                            const Transform3<T>& a2w,
                                            const Transform3<T>& b2w,
                                            const uint           method)
{
    Convex<T> const* convexA = rbA.getConvex();
    Convex<T> const* convexB = rbB.getConvex();

    Vector3<T> ptA, ptB;
    uint       nbIterGJK = 0;
    T          distance  = 0;
    distance             = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(*convexA,
                                                               *convexB,
                                                               a2w,
                                                               b2w,
                                                               rbA.getCrustThickness(),
                                                               rbB.getCrustThickness(),
                                                               ptA,
                                                               ptB,
                                                               nbIterGJK);
    return (distance);
}

// -------------------------------------------------------------------------------------------------
/** @brief RigidBody-free version: all contact properties read from ShapeData directly.
    Used by the prebuilt GPU path (BVType always OFF; BV filter ran separately).
    Rectangle-Particle: replaces virtual convexB.support with device_support_raw and
    isInside with direct params[] bounds check.
    @param sdA pre-built ShapeData for body A (shapeId-indexed)
    @param sdB pre-built ShapeData for body B (shapeId-indexed)
    @param v_b2a position of B in A-local frame
    @param q_b2a rotation of B in A-local frame
    @param contactInfo output contact information */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    auto [crustA, circRadiusA, crustB, circRadiusB, averageMass, averageRadius, contactHash]
        = getPropertiesForContact(sdA, sdB);

    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    const ConvexType typeA = sdA.type;
    const ConvexType typeB = sdB.type;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - circRadiusA - circRadiusB;
        snapshot.contactPoint
            = (circRadiusA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Rectangle-Particle Case (vtable-free via device_support_raw + params[] bounds check)
    else if(typeA == ConvexType::RECTANGLE)
    {
        Vector3<T>       ptA, ptB;
        const Vector3<T> r = v_b2a[Z] > 0 ? Vector3<T>(0, 0, -1) : Vector3<T>(0, 0, 1);
        ptB                = device_support_raw(sdB, q_b2a << r);
        transform(q_b2a, v_b2a, ptB);
        if(ptB[Z] < T(0))
        {
            ptA = Vector3<T>(ptB[X], ptB[Y], T(0));
            if(ptA[X] >= -sdA.params[0] && ptA[X] <= sdA.params[0] && ptA[Y] >= -sdA.params[1]
               && ptA[Y] <= sdA.params[1])
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case: vtable-free GJK via ShapeData
    else
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        transform(q_b2a, v_b2a, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies using pre-built ShapeData
    for vtable-free GJK support evaluation -- relative transformation (vec/quat).
    The General Case replaces virtual dispatch with device_support() over ShapeData.
    Sphere-Sphere and Rectangle-Particle paths fall through to the standard virtual variants.
    @param rbA first rigid body
    @param rbB second rigid body
    @param sdA pre-built ShapeData for rbA (slot-indexed)
    @param sdB pre-built ShapeData for rbB (slot-indexed)
    @param v_b2a position describing convex B in the A's reference frame
    @param q_b2a rotation describing convex B in the A's reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_b2a,
                                                    const Quaternion<T>& q_b2a,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2a) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_b2a,
                                         q_b2a))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        snapshot.contactVector   = v_b2a;
        snapshot.overlapDistance = norm(snapshot.contactVector) - rA - rB;
        snapshot.contactPoint    = (rA + T(.5) * snapshot.overlapDistance) * snapshot.contactVector;
        snapshot.contactVector.normalize();
    }
    // Rectangle-Particle Case (falls through to virtual supporting convex)
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        Vector3<T>       ptA, ptB;
        const Vector3<T> r = v_b2a[Z] > 0 ? Vector3<T>(0, 0, -1) : Vector3<T>(0, 0, 1);
        ptB                = convexB.support(q_b2a << r);
        transform(q_b2a, v_b2a, ptB);
        if(ptB[Z] < T(0))
        {
            ptA = Vector3<T>(ptB[X], ptB[Y], T(0));
            if(convexA.isInside(ptA))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case: vtable-free via ShapeData
    else
    {
        Vector3<T> ptA, ptB;
        uint       nbIterGJK     = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_b2a,
                                                                                   q_b2a,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        // transform(q_a2a, v_a2a, ptA);
        transform(q_b2a, v_b2a, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
    return;
}

// -------------------------------------------------------------------------------------------------
/** @brief RigidBody-free world-frame version: all contact properties from ShapeData.
    Used by the prebuilt GPU path (BVType always OFF; BV filter ran separately).
    @param sdA pre-built ShapeData for body A (shapeId-indexed)
    @param sdB pre-built ShapeData for body B (shapeId-indexed)
    @param v_a2w world-frame position of A
    @param v_b2w world-frame position of B
    @param q_a2w world-frame orientation of A
    @param q_b2w world-frame orientation of B
    @param contactInfo output contact information */
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    auto [crustA, circRadiusA, crustB, circRadiusB, averageMass, averageRadius, contactHash]
        = getPropertiesForContact(sdA, sdB);

    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    const ConvexType typeA = sdA.type;
    const ConvexType typeB = sdB.type;
    Vector3<T>       ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - circRadiusA - circRadiusB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (circRadiusA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Rectangle-Particle Case (vtable-free via device_support_raw + params[] bounds check)
    else if(typeA == ConvexType::RECTANGLE)
    {
        Vector3<T> r = q_a2w >> Vector3<T>(0, 0, 1);
        r.normalize();
        r *= copysign(T(1), r * (v_b2w - v_a2w));
        ptB = (q_b2w >> device_support_raw(sdB, q_b2w << (-r))) + v_b2w;
        if(r * (ptB - v_a2w) < T(0))
        {
            ptA                  = ((v_a2w - ptB) * r) * r + ptB;
            const Vector3<T> loc = q_a2w << (ptA - v_a2w);
            if(loc[X] >= -sdA.params[0] && loc[X] <= sdA.params[0] && loc[Y] >= -sdA.params[1]
               && loc[Y] <= sdA.params[1])
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case: vtable-free GJK via ShapeData
    else
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        transform(q_a2w, v_a2w, ptA);
        transform(q_b2w, v_b2w, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    contactInfo.setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the contact information (if any) for 2 rigid bodies using pre-built ShapeData
    for vtable-free GJK support evaluation -- world-frame (vec/quat).
    The General Case replaces virtual dispatch with device_support() over ShapeData.
    Sphere-Sphere and Rectangle-Particle paths fall through to the standard virtual variants.
    @param rbA first rigid body
    @param rbB second rigid body
    @param sdA pre-built ShapeData for rbA (slot-indexed)
    @param sdB pre-built ShapeData for rbB (slot-indexed)
    @param v_a2w position describing convex A in the world reference frame
    @param v_b2w position describing convex B in the world reference frame
    @param q_a2w rotation describing convex A in the world reference frame
    @param q_b2w rotation describing convex B in the world reference frame
    @param contactInfo output contact information */
template <typename T,
          GJKType            GJKVARIANT,
          bool               GJKACC,
          BoundingVolumeType BVType = BoundingVolumeType::OFF>
__HOSTDEVICE__ inline void closestPointsRigidBodies(const RigidBody<T>&  rbA,
                                                    const RigidBody<T>&  rbB,
                                                    const ShapeData<T>&  sdA,
                                                    const ShapeData<T>&  sdB,
                                                    const Vector3<T>&    v_a2w,
                                                    const Vector3<T>&    v_b2w,
                                                    const Quaternion<T>& q_a2w,
                                                    const Quaternion<T>& q_b2w,
                                                    ContactInfo<T>&      contactInfo)
{
    // BV pre-filter: sphere reject then OBB or OBC SAT before entering GJK
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_a2w,
                                         v_b2w,
                                         q_a2w,
                                         q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
        if(norm2(v_b2w - v_a2w) >= radiiSum * radiiSum)
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
        auto axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        const Vector3<T> bcA = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB = rbB.getConvex()->computeBoundingCylinder();
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_a2w,
                                              v_b2w,
                                              q_a2w,
                                              q_b2w))
        {
            contactInfo.setOverlapDistance(T(1));
            return;
        }
    }
    // Extract properties for contact and populate snapshot
    auto [convexA_ptr,
          crustA,
          circRadiusA,
          convexB_ptr,
          crustB,
          circRadiusB,
          averageMass,
          averageRadius,
          contactHash]
        = getPropertiesForContact(rbA, rbB);
    typename ContactInfo<T>::Snapshot snapshot;
    snapshot.averageMass     = averageMass;
    snapshot.averageRadius   = averageRadius;
    snapshot.contactHash     = contactHash;
    snapshot.overlapDistance = std::numeric_limits<T>::max();

    // Get convexes and their types
    const Convex<T>& convexA = *(convexA_ptr);
    const Convex<T>& convexB = *(convexB_ptr);
    const ConvexType typeA   = convexA.getConvexType();
    const ConvexType typeB   = convexB.getConvexType();

    Vector3<T> ptA, ptB;

    // Sphere-Sphere Case
    if(typeA == ConvexType::SPHERE && typeB == ConvexType::SPHERE)
    {
        T rA                     = rbA.getCircumscribedRadius();
        T rB                     = rbB.getCircumscribedRadius();
        ptA                      = v_a2w;
        ptB                      = v_b2w - ptA;
        snapshot.overlapDistance = norm(ptB) - rA - rB;
        {
            Vector3<T> d_hat       = ptB.normalized();
            snapshot.contactPoint  = ptA + (rA + T(.5) * snapshot.overlapDistance) * d_hat;
            snapshot.contactVector = d_hat;
        }
    }
    // Rectangle-Particle Case
    else if(typeA == ConvexType::RECTANGLE)
    {
        // ptA is the point on rectangle and ptB is the point on particle
        Vector3<T> r = q_a2w >> Vector3<T>(0, 0, 1);
        r.normalize();
        r *= copysign(T(1), r * (v_b2w - v_a2w));
        ptB = (q_b2w >> convexB.support(q_b2w << (-r))) + v_b2w;
        if(r * (ptB - v_a2w) < T(0))
        {
            ptA = ((v_a2w - ptB) * r) * r + ptB;
            if(convexA.isInside(q_a2w << (ptA - v_a2w)))
            {
                snapshot.contactPoint = T(0.5) * (ptA + ptB);
                ptB -= ptA;
                snapshot.overlapDistance = -norm(ptB);
                snapshot.contactVector   = ptB / snapshot.overlapDistance;
            }
        }
    }
    else if(typeB == ConvexType::RECTANGLE)
    {
        GAbort("General-Rectangle collision detection is not implemented yet.");
    }
    // General Case: vtable-free via ShapeData
    else
    {
        uint nbIterGJK           = 0;
        snapshot.overlapDistance = computeClosestPoints_GJK<T, GJKVARIANT, GJKACC>(sdA,
                                                                                   sdB,
                                                                                   v_a2w,
                                                                                   v_b2w,
                                                                                   q_a2w,
                                                                                   q_b2w,
                                                                                   crustA,
                                                                                   crustB,
                                                                                   ptA,
                                                                                   ptB,
                                                                                   nbIterGJK);
        if(fabs(snapshot.overlapDistance) < HIGHEPS<T>)
            snapshot.overlapDistance = -(crustA + crustB);
        else
            snapshot.overlapDistance -= crustA + crustB;
        transform(q_a2w, v_a2w, ptA);
        transform(q_b2w, v_b2w, ptB);
        snapshot.contactPoint  = T(0.5) * (ptA + ptB);
        snapshot.contactVector = (ptB - ptA).normalized();
    }

    // Set contact information
    contactInfo.setSnapshot(snapshot);
}
//@}

#endif
