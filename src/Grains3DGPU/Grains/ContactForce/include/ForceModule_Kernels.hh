// TODO: CHANGE THE FORMAT FROM HH TO CUH LATER.
#ifndef _FORCEMODULE_KERNELS_CUH_
#define _FORCEMODULE_KERNELS_CUH_

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include "BodyTag.hh"
#include "ContactForceModel.hh"
#include "ForceModuleCommon.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief GPU kernels for the ForceModule class.

    Contains the contact-force and external-force kernels.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name ForceModule GPU kernels */
//@{
// TODO: Fuse flagActivePairs_Kernel + DeviceSelect::Flagged into a single DeviceSelect::If pass
//       using an OverlapNegativeSelector predicate. The predicate approach compiles and links
//       correctly but contacts are not detected (nActive always 0).  Root cause not yet
//       identified -- likely related to ContactInfo device function visibility inside CUB's
//       internally-instantiated kernel. Revisit when upgrading CUDA/CUB version.
// NOTE: getOverlapDistance() and getSnapshot() are now inlined in ContactInfo.hh so the full
//       body is visible to CUB's template instantiation in this TU.

// -------------------------------------------------------------------------------------------------
/** @brief Predicate for CUB DeviceSelect::If that selects pairs whose overlap distance < 0.
    The counting input iterator feeds pair indices; this predicate reads the correspondng
    ContactInfo slot and returns true only for true contacts. */
template <typename T>
struct OverlapNegativeSelector
{
    const ContactInfo<T>*    contactInfo;
    __host__ __device__ bool operator()(uint pairIdx) const
    {
        return contactInfo[pairIdx].getOverlapDistance() < T(0);
    }
};

// -------------------------------------------------------------------------------------------------
/** @brief Queries the temporary storage size required by CUB DeviceSelect::If.
    @param nPairs         total number of pairs
    @param activeIdxDev   device array receiving compact indices (type deduction only)
    @param numSelectedDev device pointer receiving the count of selected items (type deduction only)
    @return required temporary storage size in bytes */
template <typename T>
INLINE size_t queryCubSelectIfTempStorageBytes(uint  nPairs,
                                               uint* activeIdxDev,
                                               uint* numSelectedDev)
{
    cub::CountingInputIterator<uint> countIter(0u);
    OverlapNegativeSelector<T>       selector{nullptr};  // nullptr safe for size query
    size_t                           bytes = 0;
    cudaErrCheck(cub::DeviceSelect::If(nullptr,
                                       bytes,
                                       countIter,
                                       activeIdxDev,
                                       numSelectedDev,
                                       static_cast<int>(nPairs),
                                       selector));
    return bytes;
}

// -------------------------------------------------------------------------------------------------
/** @brief Build a compact list of active pair indices using CUB DeviceSelect::If.
    One CUB pass: the counting iterator feeds indices 0..nPairs-1; the OverlapNegativeSelector
    predicate reads contactInfo on the device and filters in-contact pairs directly.
    @param contactInfoDev  device array of ContactInfo (length >= nPairs)
    @param nPairs          total number of pairs
    @param activeIdxDev    device array with capacity >= nPairs to receive active indices
    @param mappedCountDev  mapped-pinned device alias (getDeviceData()) to receive the count;
                           CUB writes here directly into pinned host memory -- no cudaMemcpy needed
    @param tempStorage     preallocated CUB temporary storage
    @param tempStorageBytes size of tempStorage in bytes
    @param mappedCountHost mapped-pinned host pointer (getData()) read after synchronize
    @return number of active pairs */
template <typename T>
INLINE uint buildCompactActiveIndexIf(const ContactInfo<T>* contactInfoDev,
                                      uint                  nPairs,
                                      uint*                 activeIdxDev,
                                      uint*                 mappedCountDev,
                                      void*                 tempStorage,
                                      size_t                tempStorageBytes,
                                      const uint*           mappedCountHost)
{
    cub::CountingInputIterator<uint> countIter(0u);
    OverlapNegativeSelector<T>       selector{contactInfoDev};
    cudaErrCheck(cub::DeviceSelect::If(tempStorage,
                                       tempStorageBytes,
                                       countIter,
                                       activeIdxDev,
                                       mappedCountDev,
                                       static_cast<int>(nPairs),
                                       selector));
    // CUB wrote the count directly into pinned mapped memory via the device alias;
    // synchronize so the host-side read below sees the completed value.
    cudaErrCheck(cudaDeviceSynchronize());
    return *mappedCountHost;
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the contact forces (writes to intermediate per-pair storage).
    @param CF contact force models
    @param pairList list of rigid bodies pairs
    @param contactInfo contact information
    @param activeIdx list of active pair indices
    @param position position of the components
    @param velocity kinematics of the components
    @param intermediateTorceA intermediate torce storage for particle A in each pair
    @param intermediateTorceB intermediate torce storage for particle B in each pair
    @param contactMemory view of contact memory (hash table + history data)
    @param nActive number of active pairs */
template <typename T>
__GLOBAL__ void computeContactForces_Kernel(const ContactForceModel<T>* const* CF,
                                            const uint2*                       pairList,
                                            const ContactInfo<T>*              contactInfo,
                                            const uint*                        activeIdx,
                                            const Vector3<T>*                  position,
                                            const Kinematics<T>*               velocity,
                                            Torce<T>*                          intermediateTorceA,
                                            Torce<T>*                          intermediateTorceB,
                                            ContactMemoryView<T>               contactMemory,
                                            const uint                         nActive)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;

    if(tID >= nActive)
        return;

    const uint i = activeIdx[tID];
    computeContactForces_common(CF,
                                pairList,
                                contactInfo,
                                position,
                                velocity,
                                intermediateTorceA,
                                intermediateTorceB,
                                contactMemory,
                                i);
}

// -------------------------------------------------------------------------------------------------
/** @brief Reduces per-pair intermediate torces to per-particle torces using atomics.
    @param pairList list of rigid bodies pairs
    @param activeIdx list of active pair indices
    @param intermediateTorceA intermediate torce storage for particle A in each pair
    @param intermediateTorceB intermediate torce storage for particle B in each pair
    @param torce final per-particle torce array (accumulated atomically)
    @param nActive number of active pairs */
template <typename T>
__GLOBAL__ void reduceTorces_Kernel(const uint2*    pairList,
                                    const uint*     activeIdx,
                                    const Torce<T>* intermediateTorceA,
                                    const Torce<T>* intermediateTorceB,
                                    Torce<T>*       torce,
                                    const uint      nActive)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;

    if(tID >= nActive)
        return;

    const uint i = activeIdx[tID];
    reduceTorces_common(pairList, intermediateTorceA, intermediateTorceB, torce, i);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes contact forces for ALL pairs without prior compaction.
    Each thread resets its intermediate torce slot unconditionally, then computes forces only if
    the pair is actually in contact (checked inside computeContactForces_common).
    @param CF contact force models
    @param pairList list of rigid bodies pairs
    @param contactInfo contact information
    @param position position of the components
    @param velocity kinematics of the components
    @param intermediateTorceA intermediate torce storage for particle A (indexed by pair ID)
    @param intermediateTorceB intermediate torce storage for particle B (indexed by pair ID)
    @param contactMemory view of contact memory (hash table + history data)
    @param nPairs total number of pairs */
template <typename T>
__GLOBAL__ void computeContactForces_AllPairs_Kernel(const ContactForceModel<T>* const* CF,
                                                     const uint2*                       pairList,
                                                     const ContactInfo<T>*              contactInfo,
                                                     const Vector3<T>*                  position,
                                                     const Kinematics<T>*               velocity,
                                                     Torce<T>*            intermediateTorceA,
                                                     Torce<T>*            intermediateTorceB,
                                                     ContactMemoryView<T> contactMemory,
                                                     const uint           nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;

    if(tID >= nPairs)
        return;

    // Always zero the intermediate slots so the reduce step can safely add all pairs.
    intermediateTorceA[tID].reset();
    intermediateTorceB[tID].reset();

    computeContactForces_common(CF,
                                pairList,
                                contactInfo,
                                position,
                                velocity,
                                intermediateTorceA,
                                intermediateTorceB,
                                contactMemory,
                                tID);
}

// -------------------------------------------------------------------------------------------------
/** @brief Reduces per-pair intermediate torces to per-particle torces (no-compaction variant).
    Operates over ALL pairs without an activeIdx indirection.
    @param pairList list of rigid bodies pairs
    @param intermediateTorceA intermediate torce storage for particle A in each pair
    @param intermediateTorceB intermediate torce storage for particle B in each pair
    @param torce final per-particle torce array (accumulated atomically)
    @param nPairs total number of pairs */
template <typename T>
__GLOBAL__ void reduceTorces_AllPairs_Kernel(const uint2*    pairList,
                                             const Torce<T>* intermediateTorceA,
                                             const Torce<T>* intermediateTorceB,
                                             Torce<T>*       torce,
                                             const uint      nPairs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;

    if(tID >= nPairs)
        return;

    reduceTorces_common(pairList, intermediateTorceA, intermediateTorceB, torce, tID);
}

// -------------------------------------------------------------------------------------------------
/** @brief Adds external forces such as gravity.
    @param gX x-component of the gravity vector
    @param gY y-component of the gravity vector
    @param gZ z-component of the gravity vector
    @param rigidBody array of rigid bodies for components
    @param torce array of components torces
    @param nObstacles number of obstacles
    @param nParticles number of particles */
template <typename T>
__GLOBAL__ void addExternalForces_Kernel(const T                    gX,
                                         const T                    gY,
                                         const T                    gZ,
                                         const RigidBody<T>* const* rigidBody,
                                         Torce<T>*                  torce,
                                         const uint                 nObstacles,
                                         const uint                 nParticles)
{
    uint pID = blockIdx.x * blockDim.x + threadIdx.x;

    if(pID >= nParticles)
        return;

    addExternalForces_common(Vector3<T>(gX, gY, gZ), rigidBody, torce, nObstacles + pID);
}

// -------------------------------------------------------------------------------------------------
/** @brief Accumulates forces and torques from non-master sub-bodies into their composite master,
    then resets the sub-body torces. One thread per component.
    @param torce      per-component torce (sub-bodies read, master accumulated atomically, reset)
    @param position   world-frame positions (needed to compute moment arm r = sub - master)
    @param masterSlot lookup: masterSlot[compositeIdx] = current array slot of composite master
    @param bodyTag    per-component body tag (encodes isSubBody / compositeIdx / localIdx)
    @param nComponents total number of components (obstacles + particles) */
template <typename T>
__GLOBAL__ void assembleCompositeTorces_Kernel(Torce<T>*         torce,
                                               const Vector3<T>* position,
                                               const uint*       masterSlot,
                                               const uint*       bodyTag,
                                               const uint        nComponents)
{
    uint cID = blockIdx.x * blockDim.x + threadIdx.x;
    if(cID >= nComponents)
        return;

    const uint tag = bodyTag[cID];
    if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
        return;

    const uint       mSlot = masterSlot[getCompositeIdx(tag)];
    const Vector3<T> r     = position[cID] - position[mSlot];
    const Vector3<T> f     = torce[cID].getForce();
    const Vector3<T> tau   = torce[cID].getTorque() + (r ^ f);
    torce[mSlot].addForceAtomic(f);
    torce[mSlot].addTorqueAtomic(tau);
    torce[cID].reset();
}
//@}

#endif
