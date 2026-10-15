#ifndef _COLLISIONDETECTIONMODULE_KERNELS_HH_
#define _COLLISIONDETECTIONMODULE_KERNELS_HH_

#include <cstdint>
#include <cuda_runtime.h>

#include "BodyTag.hh"
#include "ContactInfo.hh"
#include "GJK.hh"
#include "GJK_ShapeData.hh"
#include "GrainsParameters.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Transform3.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief CUDA kernel declarations for CollisionDetectionModule.

    Definitions live in CollisionDetectionModule_Kernels.cpp (compiled by nvcc).

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name CollisionDetectionModule_Kernels */
//@{

/** @brief Computes relative transformations (vec/quat) for all pairs in parallel.
    @param pairList list of pairs
    @param position world positions
    @param quaternion world quaternions
    @param relPosition output relative positions (B in A-local)
    @param relQuaternion output relative quaternions (B in A-local)
    @param nPairs number of pairs */
template <typename T>
__GLOBAL__ void computeRelativeTransformations_Kernel(const Vector3<T>*    position,
                                                      const Quaternion<T>* quaternion,
                                                      const uint2*         pairList,
                                                      Vector3<T>*          relPosition,
                                                      Quaternion<T>*       relQuaternion,
                                                      const uint           nPairs);

/** @brief Bounding volume pre-filter kernel (DEVICE path only).
    First rejects intra-composite pairs (same composite, different sub-bodies) unconditionally.
    Then applies the BV test for OBB/OBC; when BVType is OFF the BV section is compiled away,
    leaving a pure composite-membership filter reusing the same CUB compaction machinery.
    Writes a no-contact sentinel into @p contactInfo for every rejected pair.
    @param rigidBodies   Rigid body array
    @param pairList      List of pairs
    @param bodyTags      Per-component body tags (used for composite membership test)
    @param numComposites Number of composite bodies (0 = no composites, skip the check)
    @param relPosition   Per-pair relative positions (B in A-local); unused when BVType == OFF
    @param relQuaternion Per-pair relative quaternions (B in A-local); unused when BVType == OFF
    @param contactInfo   Per-pair contact info buffer; sentinel written for rejected pairs
    @param bvPassFlags   Output pass/fail flags (0 = reject, 1 = pass)
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Kernel(const RigidBody<T>* const* rigidBodies,
                                     const uint2*               pairList,
                                     const uint*                bodyTags,
                                     uint                       numComposites,
                                     const Vector3<T>*          relPosition,
                                     const Quaternion<T>*       relQuaternion,
                                     ContactInfo<T>*            contactInfo,
                                     uint8_t*                   bvPassFlags,
                                     const uint                 nPairs);

/** @brief BV pre-filter kernel using world-frame positions/quaternions directly.
    Computes relative transforms per-pair on-the-fly, eliminating the separate
    computeRelativeTransformations_Kernel pre-pass. Used when UseRelativeTransformations=false.
    @param rigidBodies   Rigid body array
    @param pairList      List of pairs
    @param bodyTags      Per-component body tags (composite membership test)
    @param numComposites Number of composite bodies (0 = no composites, skip the check)
    @param positions     World-frame particle positions
    @param quaternions   World-frame particle quaternions
    @param contactInfo   Per-pair contact info buffer; sentinel written for rejected pairs
    @param bvPassFlags   Output pass/fail flags (0 = reject, 1 = pass)
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Global_Kernel(const RigidBody<T>* const* rigidBodies,
                                            const uint2*               pairList,
                                            const uint*                bodyTags,
                                            uint                       numComposites,
                                            const Vector3<T>*          positions,
                                            const Quaternion<T>*       quaternions,
                                            ContactInfo<T>*            contactInfo,
                                            uint8_t*                   bvPassFlags,
                                            const uint                 nPairs);

/** @brief BV pre-filter kernel using pre-built BVData (circumscribed radius stored inline).
    Uses per-pair relative transforms (pre-computed by computeRelativeTransformations_Kernel).
    @param bvData        Pre-built BVData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList      List of pairs
    @param bodyTags      Per-component body tags (shapeId + composite membership)
    @param numComposites Number of composite bodies
    @param relPosition   Per-pair relative positions (B in A-local)
    @param relQuaternion Per-pair relative quaternions (B in A-local)
    @param contactInfo   Per-pair contact info buffer; sentinel written for rejected pairs
    @param bvPassFlags   Output pass/fail flags (0 = reject, 1 = pass)
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Kernel(const BVData<T>*     bvData,
                                     const uint2*         pairList,
                                     const uint*          bodyTags,
                                     uint                 numComposites,
                                     const Vector3<T>*    relPosition,
                                     const Quaternion<T>* relQuaternion,
                                     ContactInfo<T>*      contactInfo,
                                     uint8_t*             bvPassFlags,
                                     const uint           nPairs);

/** @brief BV pre-filter kernel using pre-built BVData and world-frame positions/quaternions.
    Circumscribed radii are stored in BVData; no RigidBody pointer needed.
    @param bvData        Pre-built BVData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList      List of pairs
    @param bodyTags      Per-component body tags (shapeId + composite membership)
    @param numComposites Number of composite bodies
    @param positions     World-frame positions
    @param quaternions   World-frame quaternions
    @param contactInfo   Per-pair contact info buffer; sentinel written for rejected pairs
    @param bvPassFlags   Output pass/fail flags
    @param nPairs        Number of pairs */
template <typename T, BoundingVolumeType BVType = BoundingVolumeType::OBB>
__GLOBAL__ void filterPairsBV_Global_Kernel(const BVData<T>*     bvData,
                                            const uint2*         pairList,
                                            const uint*          bodyTags,
                                            uint                 numComposites,
                                            const Vector3<T>*    positions,
                                            const Quaternion<T>* quaternions,
                                            ContactInfo<T>*      contactInfo,
                                            uint8_t*             bvPassFlags,
                                            const uint           nPairs);

/** @brief Narrow-phase GJK detection (relative vec/quat).
    When activePairIndices is non-null each thread resolves its original pair index through the
    indirection table; when null the thread ID is used directly.
    @param rigidBody         Rigid body array
    @param pairList          Full pair list
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param relPosition       Per-pair relative positions (B in A-local)
    @param relQuaternion     Per-pair relative quaternions (B in A-local)
    @param contactInfo       Output contact information (A-local frame)
    @param nPairs            Number of pairs to process (active or total) */
template <typename T, GJKType GJKVARIANT = GJKType::JOHNSON, bool GJKACC = false>
__GLOBAL__ void detectCollisionsComponents_Kernel(const RigidBody<T>* const* rigidBody,
                                                  const uint2*               pairList,
                                                  const uint*                activePairIndices,
                                                  const Vector3<T>*          relPosition,
                                                  const Quaternion<T>*       relQuaternion,
                                                  ContactInfo<T>*            contactInfo,
                                                  const uint                 nPairs);

/** @brief Narrow-phase GJK detection using absolute world-frame positions and quaternions.
    Used when relative transformations are disabled (e.g. sphere-only simulations where computing
    relative transforms is more expensive than running GJK directly in world frame).
    When @p activePairIndices is non-null each thread resolves its original pair index through the
    indirection table (BV-compacted path); when null the thread ID is used directly.
    @param rigidBody         Rigid body array
    @param pairList          Full pair list
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param position          World-frame positions
    @param quaternion        World-frame quaternions
    @param contactInfo       Output contact information
    @param nPairs            Number of pairs to process (active or total) */
template <typename T,
          GJKType            GJKVARIANT = GJKType::JOHNSON,
          bool               GJKACC     = false,
          BoundingVolumeType BVType     = BoundingVolumeType::OFF>
__GLOBAL__ void detectCollisionsComponentsGlobal_Kernel(const RigidBody<T>* const* rigidBody,
                                                        const uint2*               pairList,
                                                        const uint*          activePairIndices,
                                                        const Vector3<T>*    position,
                                                        const Quaternion<T>* quaternion,
                                                        ContactInfo<T>*      contactInfo,
                                                        const uint           nPairs);

/** @brief Transforms contact info from A-local frame to world frame for all pairs in parallel.
    @param pairList list of pairs
    @param position world positions
    @param quaternion world quaternions
    @param contactInfoLocal input contact info in A-local frame
    @param contactInfoWorld output contact info in world frame
    @param nPairs number of pairs */
template <typename T>
__GLOBAL__ void transformContactInfo_Kernel(const Vector3<T>*    position,
                                            const Quaternion<T>* quaternion,
                                            const uint2*         pairList,
                                            ContactInfo<T>*      contactInfoLocal,
                                            ContactInfo<T>*      contactInfoWorld,
                                            const uint           nPairs);
/** @brief Rebuilds the per-composite masterSlot lookup after a Morton sort.
    Only threads whose body tag identifies a sub-body with local index 0 write to the table.
    @param masterSlot  Per-composite master slot array (size = numComposites)
    @param bodyTag     Per-component body tag array (size = nComponents)
    @param nComponents Total number of components (obstacles + particles) */
__GLOBAL__ void
    rebuildMasterSlot_Kernel(uint* masterSlot, const uint* bodyTag, const uint nComponents);

/** @brief Fills the compact ShapeData table from the rigid body pointer array on the GPU.
    One thread per unique shape; repSlots[k] gives the representative component slot for shape k.
    @param shapeData      Output ShapeData array (size = nUniqueShapes)
    @param rigidBody      Rigid body pointer array (slot-indexed, size >= nComponents)
    @param repSlots       Representative slot for each unique shape (size = nUniqueShapes)
    @param nUniqueShapes  Number of unique shapes (table size) */
template <typename T>
__GLOBAL__ void fillShapeData_Kernel(ShapeData<T>*              shapeData,
                                     const RigidBody<T>* const* rigidBody,
                                     const uint*                repSlots,
                                     const uint                 nUniqueShapes);

/** @brief Fills the compact BVData table from the rigid body pointer array on the GPU.
    One thread per unique shape; repSlots[k] gives the representative component slot for shape k.
    @param bvData         Output BVData array (size = nUniqueShapes)
    @param rigidBody      Rigid body pointer array (slot-indexed, size >= nComponents)
    @param repSlots       Representative slot for each unique shape (size = nUniqueShapes)
    @param nUniqueShapes  Number of unique shapes (table size) */
template <typename T>
__GLOBAL__ void fillBVData_Kernel(BVData<T>*                 bvData,
                                  const RigidBody<T>* const* rigidBody,
                                  const uint*                repSlots,
                                  const uint                 nUniqueShapes);

/** @brief Narrow-phase GJK detection with vtable-free support evaluation via ShapeData.
    ShapeData indexed by shapeId (via bodyTags) for deduplication -- no RigidBody needed.
    When activePairIndices is non-null each thread resolves its original pair index via the
    indirection table; when null the thread ID maps directly to the pair.
    @param shapeData         Pre-built ShapeData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList          Full pair list
    @param bodyTags          Per-component body tags (encodes shapeId)
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param relPosition       Per-pair relative positions (B in A-local)
    @param relQuaternion     Per-pair relative quaternions (B in A-local)
    @param contactInfo       Output contact information (A-local frame)
    @param nPairs            Number of pairs to process (active or total) */
template <typename T, GJKType GJKVARIANT = GJKType::JOHNSON, bool GJKACC = false>
__GLOBAL__ void detectCollisionsComponents_Kernel(const ShapeData<T>*  shapeData,
                                                  const uint2*         pairList,
                                                  const uint*          bodyTags,
                                                  const uint*          activePairIndices,
                                                  const Vector3<T>*    relPosition,
                                                  const Quaternion<T>* relQuaternion,
                                                  ContactInfo<T>*      contactInfo,
                                                  const uint           nPairs);

/** @brief Narrow-phase GJK detection with vtable-free support evaluation via ShapeData,
    using absolute world-frame positions and quaternions.
    ShapeData indexed by shapeId (via bodyTags) for deduplication -- no RigidBody needed.
    When activePairIndices is non-null each thread resolves its original pair index via the
    indirection table; when null the thread ID maps directly to the pair.
    @param shapeData         Pre-built ShapeData array (shapeId-indexed, size = nUniqueShapes)
    @param pairList          Full pair list
    @param bodyTags          Per-component body tags (encodes shapeId)
    @param activePairIndices Compacted index table from CUB, or nullptr to process all pairs
    @param position          World-frame positions
    @param quaternion        World-frame quaternions
    @param contactInfo       Output contact information
    @param nPairs            Number of pairs to process (active or total) */
template <typename T, GJKType GJKVARIANT = GJKType::JOHNSON, bool GJKACC = false>
__GLOBAL__ void detectCollisionsComponentsGlobal_Kernel(const ShapeData<T>*  shapeData,
                                                        const uint2*         pairList,
                                                        const uint*          bodyTags,
                                                        const uint*          activePairIndices,
                                                        const Vector3<T>*    position,
                                                        const Quaternion<T>* quaternion,
                                                        ContactInfo<T>*      contactInfo,
                                                        const uint           nPairs);
//@}

#endif
