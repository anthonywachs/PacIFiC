#ifndef _COLLISIONDETECTIONCOMMON_HH_
#define _COLLISIONDETECTIONCOMMON_HH_

#include "CollisionDetection.hh"
#include "ContactInfo.hh"
#include "GJK_ShapeData.hh"
#include "GrainsParameters.hh"
#include "OBB.hh"
#include "OBC.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"
#include "Transform3.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief Common per-pair collision detection functions shared between CPU loops and GPU kernels.

    All functions are header-only, templated, and marked __HOSTDEVICE__ so they can be called
    from both CPU loops (in CollisionDetectionModule.cpp) and CUDA kernels.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name Collision Detection Common Functions */
//@{
/** @brief Computes relative transformations per pair (vec/quat)
    @param pairList list of rigid bodies pairs
    @param position position of the components (world frame)
    @param quaternion quaternion of the components (world frame)
    @param relativePosition output relative position of B in A-local frame
    @param relativeQuaternion output relative quaternion of B in A-local frame
    @param pairID ID of the pair */
template <typename T>
__HOSTDEVICE__ static INLINE void
    computeRelativeTransformations_common(const uint2*         pairList,
                                          const Vector3<T>*    position,
                                          const Quaternion<T>* quaternion,
                                          Vector3<T>*          relativePosition,
                                          Quaternion<T>*       relativeQuaternion,
                                          const uint           pairID)
{
    const uint2 pair           = pairList[pairID];
    const uint  idA            = pair.x;
    const uint  idB            = pair.y;
    relativePosition[pairID]   = quaternion[idA] << (position[idB] - position[idA]);
    relativeQuaternion[pairID] = inverse(quaternion[idA]) * quaternion[idB];
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes relative transformations per pair (Transform3)
    @param pairList list of rigid bodies pairs
    @param transform transforms of the components (world frame)
    @param relativeTransform output relative transform of B in A-local frame
    @param pairID ID of the pair */
template <typename T>
__HOSTDEVICE__ static INLINE void
    computeRelativeTransformations_common(const uint2*         pairList,
                                          const Transform3<T>* transform,
                                          Transform3<T>*       relativeTransform,
                                          const uint           pairID)
{
    const uint2   pair = pairList[pairID];
    const uint    idA  = pair.x;
    const uint    idB  = pair.y;
    Transform3<T> invA;
    invA.setToInverseTransform(transform[idA]);
    relativeTransform[pairID].setToTransformsComposition(transform[idB], invA);
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components (relative vec/quat).
    @param pairList list of contact pairs
    @param rigidBody rigid body array
    @param relPosition relative position of B in A-local frame
    @param relQuaternion relative quaternion of B in A-local frame
    @param contactInfo output contact information (A-local frame)
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponents_common(const uint2*               pairList,
                                      const RigidBody<T>* const* rigidBody,
                                      const Vector3<T>*          relPosition,
                                      const Quaternion<T>*       relQuaternion,
                                      ContactInfo<T>*            contactInfo,
                                      const uint                 pairID)
{
    const uint2          pair  = pairList[pairID];
    const uint           idA   = pair.x;
    const uint           idB   = pair.y;
    const RigidBody<T>&  rbA   = *(rigidBody[idA]);
    const RigidBody<T>&  rbB   = *(rigidBody[idB]);
    const Vector3<T>&    v_b2a = relPosition[pairID];
    const Quaternion<T>& q_b2a = relQuaternion[pairID];

    closestPointsRigidBodies<T, GJKVARIANT, GJKACC, BVType>(rbA,
                                                            rbB,
                                                            v_b2a,
                                                            q_b2a,
                                                            contactInfo[pairID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components using global coordinates (vec/quat)
    @param pairList list of contact pairs
    @param rigidBody rigid body array
    @param position global position of the components
    @param quaternion global quaternion of the components
    @param contactInfo output contact information
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponentsGlobal_common(const uint2*               pairList,
                                            const RigidBody<T>* const* rigidBody,
                                            const Vector3<T>*          position,
                                            const Quaternion<T>*       quaternion,
                                            ContactInfo<T>*            contactInfo,
                                            const uint                 pairID)
{
    const uint2          pair  = pairList[pairID];
    const uint           idA   = pair.x;
    const uint           idB   = pair.y;
    const RigidBody<T>&  rbA   = *(rigidBody[idA]);
    const RigidBody<T>&  rbB   = *(rigidBody[idB]);
    const Vector3<T>&    v_a2w = position[idA];
    const Vector3<T>&    v_b2w = position[idB];
    const Quaternion<T>& q_a2w = quaternion[idA];
    const Quaternion<T>& q_b2w = quaternion[idB];
    closestPointsRigidBodies<T, GJKVARIANT, GJKACC, BVType>(rbA,
                                                            rbB,
                                                            v_a2w,
                                                            v_b2w,
                                                            q_a2w,
                                                            q_b2w,
                                                            contactInfo[pairID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components (relative Transform3)
    @param pairList list of contact pairs
    @param rigidBody rigid body array
    @param relTransform relative transform of B in A-local frame
    @param contactInfo output contact information (A-local frame)
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponents_common(const uint2*               pairList,
                                      const RigidBody<T>* const* rigidBody,
                                      const Transform3<T>*       relTransform,
                                      ContactInfo<T>*            contactInfo,
                                      const uint                 pairID)
{
    const uint2          pair = pairList[pairID];
    const uint           idA  = pair.x;
    const uint           idB  = pair.y;
    const RigidBody<T>&  rbA  = *(rigidBody[idA]);
    const RigidBody<T>&  rbB  = *(rigidBody[idB]);
    const Transform3<T>& b2a  = relTransform[pairID];
    closestPointsRigidBodies<T, GJKVARIANT, GJKACC, BVType>(rbA, rbB, b2a, contactInfo[pairID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components using global coordinates (Transform3)
    @param pairList list of contact pairs
    @param rigidBody rigid body array
    @param transform global transforms of the components
    @param contactInfo output contact information
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponentsGlobal_common(const uint2*               pairList,
                                            const RigidBody<T>* const* rigidBody,
                                            const Transform3<T>*       transform,
                                            ContactInfo<T>*            contactInfo,
                                            const uint                 pairID)
{
    const uint2          pair = pairList[pairID];
    const uint           idA  = pair.x;
    const uint           idB  = pair.y;
    const RigidBody<T>&  rbA  = *(rigidBody[idA]);
    const RigidBody<T>&  rbB  = *(rigidBody[idB]);
    const Transform3<T>& a2w  = transform[idA];
    const Transform3<T>& b2w  = transform[idB];
    closestPointsRigidBodies<T, GJKVARIANT, GJKACC, BVType>(rbA,
                                                            rbB,
                                                            a2w,
                                                            b2w,
                                                            contactInfo[pairID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Transforms contact information from A-local frame to world frame.
    @param pairList list of contact pairs
    @param position world positions of components
    @param quaternion world orientations of components
    @param contactInfoLocal contact info in A-local frame (input)
    @param contactInfoWorld contact info in world frame (output, only for active pairs)
    @param pairID ID of the pair */
template <typename T>
__HOSTDEVICE__ static INLINE void transformContactInfo_common(const uint2*         pairList,
                                                              const Vector3<T>*    position,
                                                              const Quaternion<T>* quaternion,
                                                              ContactInfo<T>*      contactInfoLocal,
                                                              ContactInfo<T>*      contactInfoWorld,
                                                              const uint           pairID)
{
    uint                              idA      = pairList[pairID].x;
    const Quaternion<T>&              qA       = quaternion[idA];
    const ContactInfo<T>              ciL      = contactInfoLocal[pairID];
    typename ContactInfo<T>::Snapshot snapshot = ciL.getSnapshot();
    snapshot.contactPoint                      = (qA >> snapshot.contactPoint) + position[idA];
    snapshot.contactVector                     = qA >> snapshot.contactVector;
    contactInfoWorld[pairID].setSnapshot(snapshot);
}

// -------------------------------------------------------------------------------------------------
/** @brief BV-only test for a single pair (bounding sphere + optional OBB SAT / OBC test).
    Returns true if the pair passes all BV filters and should proceed to narrow-phase GJK.
    For rejected pairs the caller is responsible for writing a no-contact sentinel to contactInfo.
    @param rbA rigid body A
    @param rbB rigid body B
    @param v_b2a relative position of B in A-local frame
    @param q_b2a relative quaternion of B w.r.t. A */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__HOSTDEVICE__ static INLINE bool filterPairBV_common(const RigidBody<T>&  rbA,
                                                      const RigidBody<T>&  rbB,
                                                      const Vector3<T>&    v_b2a,
                                                      const Quaternion<T>& q_b2a)
{
    const T radiiSum = rbA.getCircumscribedRadius() + rbB.getCircumscribedRadius();
    if(norm2(v_b2a) >= radiiSum * radiiSum)
        return false;
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        if(!intersectOrientedBoundingBox(rbA.getConvex()->computeBoundingBox(),
                                         rbB.getConvex()->computeBoundingBox(),
                                         v_b2a,
                                         q_b2a))
            return false;
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        // bc = [radius, halfHeight, axisIndex(0=X,1=Y,2=Z)]
        const Vector3<T> bcA           = rbA.getConvex()->computeBoundingCylinder();
        const Vector3<T> bcB           = rbB.getConvex()->computeBoundingCylinder();
        auto             axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_b2a,
                                              q_b2a))
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief BV-only test for a single pair using pre-built BVData (vtable-free bounding volume
    queries). Functionally identical to filterPairBV_common but replaces
    getConvex()->computeBoundingBox() / computeBoundingCylinder() with stored BVData fields.
    Circumscribed radii are read from BVData instead of RigidBody.
    @param bvA  pre-built BVData for body A (shapeId-indexed)
    @param bvB  pre-built BVData for body B (shapeId-indexed)
    @param v_b2a relative position of B in A-local frame
    @param q_b2a relative quaternion of B w.r.t. A */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__HOSTDEVICE__ static INLINE bool filterPairBV_common(const BVData<T>&     bvA,
                                                      const BVData<T>&     bvB,
                                                      const Vector3<T>&    v_b2a,
                                                      const Quaternion<T>& q_b2a)
{
    const T radiiSum = bvA.circumscribedRadius + bvB.circumscribedRadius;
    if(norm2(v_b2a) >= radiiSum * radiiSum)
        return false;
    if constexpr(BVType == BoundingVolumeType::OBB)
    {
        if(!intersectOrientedBoundingBox(bvA.boundingBox, bvB.boundingBox, v_b2a, q_b2a))
            return false;
    }
    if constexpr(BVType == BoundingVolumeType::OBC)
    {
        // bc = {radius, halfHeight, axisIndex(0=X,1=Y,2=Z)}
        const Vector3<T>& bcA           = bvA.boundingCylinder;
        const Vector3<T>& bcB           = bvB.boundingCylinder;
        auto              axisFromIndex = [](T idx) -> Vector3<T> {
            if(idx == T(0))
                return Vector3<T>(T(1), T(0), T(0));
            if(idx == T(1))
                return Vector3<T>(T(0), T(1), T(0));
            return Vector3<T>(T(0), T(0), T(1));
        };
        if(!intersectOrientedBoundingCylinder(bcA[X],
                                              bcA[Y],
                                              axisFromIndex(bcA[Z]),
                                              bcB[X],
                                              bcB[Y],
                                              axisFromIndex(bcB[Z]),
                                              v_b2a,
                                              q_b2a))
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components using pre-built ShapeData (vtable-free GJK).
    ShapeData is indexed by shapeId (from bodyTags) for deduplication across identical shapes.
    @param pairList list of contact pairs (slot indices)
    @param shapeData pre-built ShapeData array (shapeId-indexed, size = nUniqueShapes)
    @param bodyTags body-tag array (slot-indexed); encodes shapeId via getShapeId()
    @param relPosition relative position of B in A-local frame (pre-cached)
    @param relQuaternion relative quaternion of B w.r.t. A (pre-cached)
    @param contactInfo output contact information (A-local frame)
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponents_common(const uint2*         pairList,
                                      const ShapeData<T>*  shapeData,
                                      const uint*          bodyTags,
                                      const Vector3<T>*    relPosition,
                                      const Quaternion<T>* relQuaternion,
                                      ContactInfo<T>*      contactInfo,
                                      const uint           pairID)
{
    const uint2          pair  = pairList[pairID];
    const ShapeData<T>&  sdA   = shapeData[getShapeId(bodyTags[pair.x])];
    const ShapeData<T>&  sdB   = shapeData[getShapeId(bodyTags[pair.y])];
    const Vector3<T>&    v_b2a = relPosition[pairID];
    const Quaternion<T>& q_b2a = relQuaternion[pairID];

    closestPointsRigidBodies<T, GJKVARIANT, GJKACC>(sdA, sdB, v_b2a, q_b2a, contactInfo[pairID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Detects collisions between components using pre-built ShapeData and global coords.
    ShapeData is indexed by shapeId (from bodyTags) for deduplication across identical shapes.
    @param pairList list of contact pairs
    @param shapeData pre-built ShapeData array (shapeId-indexed, size = nUniqueShapes)
    @param bodyTags body-tag array (slot-indexed); encodes shapeId via getShapeId()
    @param position global position of the components
    @param quaternion global quaternion of the components
    @param contactInfo output contact information
    @param pairID ID of the pair */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__HOSTDEVICE__ static INLINE void
    detectCollisionsComponentsGlobal_common(const uint2*         pairList,
                                            const ShapeData<T>*  shapeData,
                                            const uint*          bodyTags,
                                            const Vector3<T>*    position,
                                            const Quaternion<T>* quaternion,
                                            ContactInfo<T>*      contactInfo,
                                            const uint           pairID)
{
    const uint2          pair  = pairList[pairID];
    const uint           idA   = pair.x;
    const uint           idB   = pair.y;
    const ShapeData<T>&  sdA   = shapeData[getShapeId(bodyTags[idA])];
    const ShapeData<T>&  sdB   = shapeData[getShapeId(bodyTags[idB])];
    const Vector3<T>&    v_a2w = position[idA];
    const Vector3<T>&    v_b2w = position[idB];
    const Quaternion<T>& q_a2w = quaternion[idA];
    const Quaternion<T>& q_b2w = quaternion[idB];

    closestPointsRigidBodies<T, GJKVARIANT, GJKACC>(sdA,
                                                    sdB,
                                                    v_a2w,
                                                    v_b2w,
                                                    q_a2w,
                                                    q_b2w,
                                                    contactInfo[pairID]);
}
//@}

#endif
