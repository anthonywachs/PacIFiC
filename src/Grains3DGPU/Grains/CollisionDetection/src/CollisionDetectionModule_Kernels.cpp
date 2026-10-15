#include "CollisionDetectionModule_Kernels.hh"
#include "CollisionDetectionCommon.hh"
#include <cstdint>

// -------------------------------------------------------------------------------------------------
// Computes per-pair relative position / quaternion (B in A-local frame)
template <typename T>
__GLOBAL__ void computeRelativeTransformations_Kernel(const Vector3<T>*    position,
                                                      const Quaternion<T>* quaternion,
                                                      const uint2*         pairList,
                                                      Vector3<T>*          relPosition,
                                                      Quaternion<T>*       relQuaternion,
                                                      const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    computeRelativeTransformations_common(pairList,
                                          position,
                                          quaternion,
                                          relPosition,
                                          relQuaternion,
                                          tID);
}

// -------------------------------------------------------------------------------------------------
// BV pre-filter kernel.
// Rejects intra-composite pairs unconditionally, then applies the BV test.
// When BVType is OFF the BV section is compiled away, leaving a pure composite-membership test.
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Kernel(const RigidBody<T>* const* rigidBodies,
                                     const uint2*               pairList,
                                     const uint*                bodyTags,
                                     uint                       numComposites,
                                     const Vector3<T>*          relPosition,
                                     const Quaternion<T>*       relQuaternion,
                                     ContactInfo<T>*            contactInfo,
                                     uint8_t*                   bvPassFlags,
                                     const uint                 nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    // Reject intra-composite pairs (same composite, different sub-bodies)
    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        // No BV check: all non-composite pairs pass
        bvPassFlags[tID] = 1;
    }
    else
    {
        const RigidBody<T>& rbA = *(rigidBodies[pair.x]);
        const RigidBody<T>& rbB = *(rigidBodies[pair.y]);
        if(filterPairBV_common<T, BVType>(rbA, rbB, relPosition[tID], relQuaternion[tID]))
        {
            bvPassFlags[tID] = 1;
        }
        else
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
        }
    }
}

// -------------------------------------------------------------------------------------------------
// BV pre-filter using world-frame positions/quaternions; computes relative transforms per-pair
// on-the-fly to avoid a separate computeRelativeTransformations_Kernel pre-pass.
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Global_Kernel(const RigidBody<T>* const* rigidBodies,
                                            const uint2*               pairList,
                                            const uint*                bodyTags,
                                            uint                       numComposites,
                                            const Vector3<T>*          positions,
                                            const Quaternion<T>*       quaternions,
                                            ContactInfo<T>*            contactInfo,
                                            uint8_t*                   bvPassFlags,
                                            const uint                 nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    // Reject intra-composite pairs (same composite, different sub-bodies)
    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        bvPassFlags[tID] = 1;
    }
    else
    {
        const RigidBody<T>& rbA   = *(rigidBodies[pair.x]);
        const RigidBody<T>& rbB   = *(rigidBodies[pair.y]);
        const Vector3<T>    v_b2a = quaternions[pair.x] << (positions[pair.y] - positions[pair.x]);
        const Quaternion<T> q_b2a = inverse(quaternions[pair.x]) * quaternions[pair.y];
        if(filterPairBV_common<T, BVType>(rbA, rbB, v_b2a, q_b2a))
        {
            bvPassFlags[tID] = 1;
        }
        else
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection (relative vec/quat).
// When activePairIndices is non-null each thread resolves its pair index through the indirection
// table (BV-compacted path). When null the thread ID is used directly (BV-off path, all pairs).
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__GLOBAL__ void detectCollisionsComponents_Kernel(const RigidBody<T>* const* rigidBody,
                                                  const uint2*               pairList,
                                                  const uint*                activePairIndices,
                                                  const Vector3<T>*          relPosition,
                                                  const Quaternion<T>*       relQuaternion,
                                                  ContactInfo<T>*            contactInfo,
                                                  const uint                 nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    detectCollisionsComponents_common<T, GJKVARIANT, GJKACC, BoundingVolumeType::OFF>(pairList,
                                                                                      rigidBody,
                                                                                      relPosition,
                                                                                      relQuaternion,
                                                                                      contactInfo,
                                                                                      pairIdx);
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection using absolute world-frame positions and quaternions.
// When activePairIndices is non-null each thread resolves its pair index through the indirection
// table (BV-compacted path). When null the thread ID is used directly (BV-off path, all pairs).
template <typename T, GJKType GJKVARIANT, bool GJKACC, BoundingVolumeType BVType>
__GLOBAL__ void detectCollisionsComponentsGlobal_Kernel(const RigidBody<T>* const* rigidBody,
                                                        const uint2*               pairList,
                                                        const uint*          activePairIndices,
                                                        const Vector3<T>*    position,
                                                        const Quaternion<T>* quaternion,
                                                        ContactInfo<T>*      contactInfo,
                                                        const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    detectCollisionsComponentsGlobal_common<T, GJKVARIANT, GJKACC, BVType>(pairList,
                                                                           rigidBody,
                                                                           position,
                                                                           quaternion,
                                                                           contactInfo,
                                                                           pairIdx);
}

// -------------------------------------------------------------------------------------------------
// Rebuilds masterSlot lookup after a Morton sort
__GLOBAL__ void
    rebuildMasterSlot_Kernel(uint* masterSlot, const uint* bodyTag, const uint nComponents)
{
    uint cID = blockIdx.x * blockDim.x + threadIdx.x;
    if(cID >= nComponents)
        return;

    const uint tag = bodyTag[cID];
    if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
        masterSlot[getCompositeIdx(tag)] = cID;
}

// -------------------------------------------------------------------------------------------------
// Transforms contact info from A-local frame to world frame
template <typename T>
__GLOBAL__ void transformContactInfo_Kernel(const Vector3<T>*    position,
                                            const Quaternion<T>* quaternion,
                                            const uint2*         pairList,
                                            ContactInfo<T>*      contactInfoLocal,
                                            ContactInfo<T>*      contactInfoWorld,
                                            const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    transformContactInfo_common(pairList,
                                position,
                                quaternion,
                                contactInfoLocal,
                                contactInfoWorld,
                                tID);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations
#define X(T)                                                                                \
    template __GLOBAL__ void computeRelativeTransformations_Kernel<T>(const Vector3<T>*,    \
                                                                      const Quaternion<T>*, \
                                                                      const uint2*,         \
                                                                      Vector3<T>*,          \
                                                                      Quaternion<T>*,       \
                                                                      const uint);          \
    template __GLOBAL__ void transformContactInfo_Kernel<T>(const Vector3<T>*,              \
                                                            const Quaternion<T>*,           \
                                                            const uint2*,                   \
                                                            ContactInfo<T>*,                \
                                                            ContactInfo<T>*,                \
                                                            const uint);
X(float)
X(double)
#undef X

#define X(T, GJK, ACC)                                                       \
    template __GLOBAL__ void detectCollisionsComponents_Kernel<T, GJK, ACC>( \
        const RigidBody<T>* const*,                                          \
        const uint2*,                                                        \
        const uint*,                                                         \
        const Vector3<T>*,                                                   \
        const Quaternion<T>*,                                                \
        ContactInfo<T>*,                                                     \
        const uint);
X(float, GJKType::JOHNSON, false)
X(double, GJKType::JOHNSON, false)
X(float, GJKType::JOHNSON, true)
X(double, GJKType::JOHNSON, true)
X(float, GJKType::SIGNEDVOLUME, false)
X(double, GJKType::SIGNEDVOLUME, false)
X(float, GJKType::SIGNEDVOLUME, true)
X(double, GJKType::SIGNEDVOLUME, true)
#undef X

#define X(T, GJK, ACC, BV)                                                             \
    template __GLOBAL__ void detectCollisionsComponentsGlobal_Kernel<T, GJK, ACC, BV>( \
        const RigidBody<T>* const*,                                                    \
        const uint2*,                                                                  \
        const uint*,                                                                   \
        const Vector3<T>*,                                                             \
        const Quaternion<T>*,                                                          \
        ContactInfo<T>*,                                                               \
        const uint);
X(float, GJKType::JOHNSON, false, BoundingVolumeType::OFF)
X(double, GJKType::JOHNSON, false, BoundingVolumeType::OFF)
X(float, GJKType::JOHNSON, true, BoundingVolumeType::OFF)
X(double, GJKType::JOHNSON, true, BoundingVolumeType::OFF)
X(float, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OFF)
X(double, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OFF)
X(float, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OFF)
X(double, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OFF)
X(float, GJKType::JOHNSON, false, BoundingVolumeType::OBB)
X(double, GJKType::JOHNSON, false, BoundingVolumeType::OBB)
X(float, GJKType::JOHNSON, true, BoundingVolumeType::OBB)
X(double, GJKType::JOHNSON, true, BoundingVolumeType::OBB)
X(float, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OBB)
X(double, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OBB)
X(float, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OBB)
X(double, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OBB)
X(float, GJKType::JOHNSON, false, BoundingVolumeType::OBC)
X(double, GJKType::JOHNSON, false, BoundingVolumeType::OBC)
X(float, GJKType::JOHNSON, true, BoundingVolumeType::OBC)
X(double, GJKType::JOHNSON, true, BoundingVolumeType::OBC)
X(float, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OBC)
X(double, GJKType::SIGNEDVOLUME, false, BoundingVolumeType::OBC)
X(float, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OBC)
X(double, GJKType::SIGNEDVOLUME, true, BoundingVolumeType::OBC)
#undef X

#define X(T, BV)                                                                     \
    template __GLOBAL__ void filterPairsBV_Kernel<T, BV>(const RigidBody<T>* const*, \
                                                         const uint2*,               \
                                                         const uint*,                \
                                                         uint,                       \
                                                         const Vector3<T>*,          \
                                                         const Quaternion<T>*,       \
                                                         ContactInfo<T>*,            \
                                                         uint8_t*,                   \
                                                         const uint);
X(float, BoundingVolumeType::OBB)
X(double, BoundingVolumeType::OBB)
X(float, BoundingVolumeType::OBC)
X(double, BoundingVolumeType::OBC)
X(float, BoundingVolumeType::OFF)
X(double, BoundingVolumeType::OFF)
#undef X

#define X(T, BV)                                                                            \
    template __GLOBAL__ void filterPairsBV_Global_Kernel<T, BV>(const RigidBody<T>* const*, \
                                                                const uint2*,               \
                                                                const uint*,                \
                                                                uint,                       \
                                                                const Vector3<T>*,          \
                                                                const Quaternion<T>*,       \
                                                                ContactInfo<T>*,            \
                                                                uint8_t*,                   \
                                                                const uint);
X(float, BoundingVolumeType::OBB)
X(double, BoundingVolumeType::OBB)
X(float, BoundingVolumeType::OBC)
X(double, BoundingVolumeType::OBC)
X(float, BoundingVolumeType::OFF)
X(double, BoundingVolumeType::OFF)
#undef X

// -------------------------------------------------------------------------------------------------
// BV pre-filter using pre-built BVData (vtable-free bounding-volume queries, relative frame).
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Kernel(const BVData<T>*     bvData,
                                     const uint2*         pairList,
                                     const uint*          bodyTags,
                                     uint                 numComposites,
                                     const Vector3<T>*    relPosition,
                                     const Quaternion<T>* relQuaternion,
                                     ContactInfo<T>*      contactInfo,
                                     uint8_t*             bvPassFlags,
                                     const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        bvPassFlags[tID] = 1;
    }
    else
    {
        const BVData<T>& bvA = bvData[getShapeId(bodyTags[pair.x])];
        const BVData<T>& bvB = bvData[getShapeId(bodyTags[pair.y])];
        if(filterPairBV_common<T, BVType>(bvA, bvB, relPosition[tID], relQuaternion[tID]))
        {
            bvPassFlags[tID] = 1;
        }
        else
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
        }
    }
}

// -------------------------------------------------------------------------------------------------
// BV pre-filter using pre-built BVData and world-frame transforms (vtable-free).
template <typename T, BoundingVolumeType BVType>
__GLOBAL__ void filterPairsBV_Global_Kernel(const BVData<T>*     bvData,
                                            const uint2*         pairList,
                                            const uint*          bodyTags,
                                            uint                 numComposites,
                                            const Vector3<T>*    positions,
                                            const Quaternion<T>* quaternions,
                                            ContactInfo<T>*      contactInfo,
                                            uint8_t*             bvPassFlags,
                                            const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint2 pair = pairList[tID];

    if(numComposites > 0)
    {
        const uint tagA = bodyTags[pair.x];
        const uint tagB = bodyTags[pair.y];
        if(isSubBody(tagA) && isSubBody(tagB) && getCompositeIdx(tagA) == getCompositeIdx(tagB))
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
            return;
        }
    }

    if constexpr(BVType == BoundingVolumeType::OFF)
    {
        bvPassFlags[tID] = 1;
    }
    else
    {
        const BVData<T>&    bvA   = bvData[getShapeId(bodyTags[pair.x])];
        const BVData<T>&    bvB   = bvData[getShapeId(bodyTags[pair.y])];
        const Vector3<T>    v_b2a = quaternions[pair.x] << (positions[pair.y] - positions[pair.x]);
        const Quaternion<T> q_b2a = inverse(quaternions[pair.x]) * quaternions[pair.y];
        if(filterPairBV_common<T, BVType>(bvA, bvB, v_b2a, q_b2a))
        {
            bvPassFlags[tID] = 1;
        }
        else
        {
            bvPassFlags[tID] = 0;
            contactInfo[tID].setOverlapDistance(T(1));
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations for Prebuilt BV filter kernels
#define X(T, BV)                                                               \
    template __GLOBAL__ void filterPairsBV_Kernel<T, BV>(const BVData<T>*,     \
                                                         const uint2*,         \
                                                         const uint*,          \
                                                         uint,                 \
                                                         const Vector3<T>*,    \
                                                         const Quaternion<T>*, \
                                                         ContactInfo<T>*,      \
                                                         uint8_t*,             \
                                                         const uint);
X(float, BoundingVolumeType::OBB)
X(double, BoundingVolumeType::OBB)
X(float, BoundingVolumeType::OBC)
X(double, BoundingVolumeType::OBC)
X(float, BoundingVolumeType::OFF)
X(double, BoundingVolumeType::OFF)
#undef X

#define X(T, BV)                                                                      \
    template __GLOBAL__ void filterPairsBV_Global_Kernel<T, BV>(const BVData<T>*,     \
                                                                const uint2*,         \
                                                                const uint*,          \
                                                                uint,                 \
                                                                const Vector3<T>*,    \
                                                                const Quaternion<T>*, \
                                                                ContactInfo<T>*,      \
                                                                uint8_t*,             \
                                                                const uint);
X(float, BoundingVolumeType::OBB)
X(double, BoundingVolumeType::OBB)
X(float, BoundingVolumeType::OBC)
X(double, BoundingVolumeType::OBC)
X(float, BoundingVolumeType::OFF)
X(double, BoundingVolumeType::OFF)
#undef X

// -------------------------------------------------------------------------------------------------
// Fills the compact ShapeData table; one thread per unique shape.
// Thread k reads rigidBody[repSlots[k]] and fills shapeData[k].
template <typename T>
__GLOBAL__ void fillShapeData_Kernel(ShapeData<T>*              shapeData,
                                     const RigidBody<T>* const* rigidBody,
                                     const uint*                repSlots,
                                     const uint                 nUniqueShapes)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nUniqueShapes)
        return;
    fillShapeData(shapeData[tID], rigidBody[repSlots[tID]]);
}

// -------------------------------------------------------------------------------------------------
// Fills the compact BVData table; one thread per unique shape.
// Thread k reads rigidBody[repSlots[k]] and fills bvData[k].
template <typename T>
__GLOBAL__ void fillBVData_Kernel(BVData<T>*                 bvData,
                                  const RigidBody<T>* const* rigidBody,
                                  const uint*                repSlots,
                                  const uint                 nUniqueShapes)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nUniqueShapes)
        return;
    fillBVData(bvData[tID], rigidBody[repSlots[tID]]);
}

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection with vtable-free support evaluation via ShapeData.
// ShapeData indexed by shapeId (via bodyTags); no RigidBody pointer needed.
// When activePairIndices is non-null each thread resolves its pair index through the indirection
// table (BV-compacted path). When null the thread ID is used directly (BV-off path, all pairs).
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__GLOBAL__ void detectCollisionsComponents_Kernel(const ShapeData<T>*  shapeData,
                                                  const uint2*         pairList,
                                                  const uint*          bodyTags,
                                                  const uint*          activePairIndices,
                                                  const Vector3<T>*    relPosition,
                                                  const Quaternion<T>* relQuaternion,
                                                  ContactInfo<T>*      contactInfo,
                                                  const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    detectCollisionsComponents_common<T, GJKVARIANT, GJKACC, BoundingVolumeType::OFF>(pairList,
                                                                                      shapeData,
                                                                                      bodyTags,
                                                                                      relPosition,
                                                                                      relQuaternion,
                                                                                      contactInfo,
                                                                                      pairIdx);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations for ShapeData kernels
#define X(T)                                                                     \
    template __GLOBAL__ void fillShapeData_Kernel<T>(ShapeData<T>*,              \
                                                     const RigidBody<T>* const*, \
                                                     const uint*,                \
                                                     const uint);
X(float)
X(double)
#undef X

#define X(T)                                                                  \
    template __GLOBAL__ void fillBVData_Kernel<T>(BVData<T>*,                 \
                                                  const RigidBody<T>* const*, \
                                                  const uint*,                \
                                                  const uint);
X(float)
X(double)
#undef X

#define X(T, GJK, ACC)                                                                            \
    template __GLOBAL__ void detectCollisionsComponents_Kernel<T, GJK, ACC>(const ShapeData<T>*,  \
                                                                            const uint2*,         \
                                                                            const uint*,          \
                                                                            const uint*,          \
                                                                            const Vector3<T>*,    \
                                                                            const Quaternion<T>*, \
                                                                            ContactInfo<T>*,      \
                                                                            const uint);
X(float, GJKType::JOHNSON, false)
X(double, GJKType::JOHNSON, false)
X(float, GJKType::JOHNSON, true)
X(double, GJKType::JOHNSON, true)
X(float, GJKType::SIGNEDVOLUME, false)
X(double, GJKType::SIGNEDVOLUME, false)
X(float, GJKType::SIGNEDVOLUME, true)
X(double, GJKType::SIGNEDVOLUME, true)
#undef X

// -------------------------------------------------------------------------------------------------
// Narrow-phase GJK detection with vtable-free support evaluation via ShapeData (world frame).
// ShapeData indexed by shapeId (via bodyTags); no RigidBody pointer needed.
// When activePairIndices is non-null each thread resolves its pair index through the indirection
// table (BV-compacted path). When null the thread ID is used directly (BV-off path, all pairs).
template <typename T, GJKType GJKVARIANT, bool GJKACC>
__GLOBAL__ void detectCollisionsComponentsGlobal_Kernel(const ShapeData<T>*  shapeData,
                                                        const uint2*         pairList,
                                                        const uint*          bodyTags,
                                                        const uint*          activePairIndices,
                                                        const Vector3<T>*    position,
                                                        const Quaternion<T>* quaternion,
                                                        ContactInfo<T>*      contactInfo,
                                                        const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= nPairs)
        return;

    const uint pairIdx = (activePairIndices != nullptr) ? activePairIndices[tID] : tID;
    detectCollisionsComponentsGlobal_common<T, GJKVARIANT, GJKACC, BoundingVolumeType::OFF>(
        pairList,
        shapeData,
        bodyTags,
        position,
        quaternion,
        contactInfo,
        pairIdx);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations for Global PrebuiltKernel
#define X(T, GJK, ACC)                                                             \
    template __GLOBAL__ void detectCollisionsComponentsGlobal_Kernel<T, GJK, ACC>( \
        const ShapeData<T>*,                                                       \
        const uint2*,                                                              \
        const uint*,                                                               \
        const uint*,                                                               \
        const Vector3<T>*,                                                         \
        const Quaternion<T>*,                                                      \
        ContactInfo<T>*,                                                           \
        const uint);
X(float, GJKType::JOHNSON, false)
X(double, GJKType::JOHNSON, false)
X(float, GJKType::JOHNSON, true)
X(double, GJKType::JOHNSON, true)
X(float, GJKType::SIGNEDVOLUME, false)
X(double, GJKType::SIGNEDVOLUME, false)
X(float, GJKType::SIGNEDVOLUME, true)
X(double, GJKType::SIGNEDVOLUME, true)
#undef X
