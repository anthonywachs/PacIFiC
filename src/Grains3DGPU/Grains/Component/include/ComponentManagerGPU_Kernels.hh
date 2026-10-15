// TODO: CHANGE THE FORMAT FROM HH TO CUH LATER.
#ifndef _COMPONENTMANAGERGPU_KERNLES_CUH_
#define _COMPONENTMANAGERGPU_KERNLES_CUH_

#include <cuda_runtime.h>

#include "ComponentManager.hh"
#include "ComponentManagerCommon.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"
#include "TimeIntegrator.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief GPU kernels for the ComponentManagerGPU class (particle motion).

    Contact-force kernels have been moved to ForceModule_Kernels.hh.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name ComponentManagerGPU_Kernels : External methods */
//@{
/** @brief Updates the position and velocities of particles (KDK Leapfrog Step 1).
    Non-master sub-bodies are skipped; their world transforms are updated after this call by
    updateSubBodyPositions_Kernel.
    @param TI time integrator scheme
    @param rigidBody array of rigid bodies for components
    @param position the position of the component
    @param quaternion array of components quaternions
    @param velocity array of components velocities
    @param torce array of components torces
    @param bodyTag per-component body tag (isSubBody|compositeIdx|localIdx)
    @param nObstacles number of obstacles
    @param nParticles number of particles */
template <typename T>
__GLOBAL__ void moveParticles_Kernel(const TimeIntegrator<T>* const* TI,
                                     const RigidBody<T>* const*      rigidBody,
                                     Vector3<T>*                     position,
                                     Quaternion<T>*                  quaternion,
                                     Kinematics<T>*                  velocity,
                                     Torce<T>*                       torce,
                                     const uint*                     bodyTag,
                                     const uint                      nObstacles,
                                     const uint                      nParticles)
{
    uint pID = blockIdx.x * blockDim.x + threadIdx.x;

    if(pID >= nParticles)
        return;

    const uint cID = nObstacles + pID;
    // Skip non-master sub-bodies -- their positions are set by updateSubBodyPositions_Kernel
    if(isSubBody(bodyTag[cID]) && getSubBodyLocalIdx(bodyTag[cID]) != 0u)
        return;

    moveParticles_common(TI, rigidBody, position, quaternion, velocity, torce, cID);
}

// -------------------------------------------------------------------------------------------------
/** @brief Performs the second velocity half-kick for split-step schemes (KDK Leapfrog Step 3).
    For single-pass schemes the underlying AdvanceVelocity is a no-op.
    @param TI time integrator scheme
    @param rigidBody array of rigid bodies for components
    @param quaternion array of components quaternions
    @param velocity array of components velocities
    @param torce array of components torces (read but NOT reset)
    @param bodyTag per-component body tag (isSubBody|compositeIdx|localIdx)
    @param nObstacles number of obstacles
    @param nParticles number of particles */
template <typename T>
__GLOBAL__ void advanceVelocity_Kernel(const TimeIntegrator<T>* const* TI,
                                       const RigidBody<T>* const*      rigidBody,
                                       const Quaternion<T>*            quaternion,
                                       Kinematics<T>*                  velocity,
                                       const Torce<T>*                 torce,
                                       const uint*                     bodyTag,
                                       const uint                      nObstacles,
                                       const uint                      nParticles)
{
    uint pID = blockIdx.x * blockDim.x + threadIdx.x;

    if(pID >= nParticles)
        return;

    const uint cID = nObstacles + pID;
    // Skip non-master sub-bodies -- they are slaved to their composite master
    if(isSubBody(bodyTag[cID]) && getSubBodyLocalIdx(bodyTag[cID]) != 0u)
        return;

    advanceVelocity_common(TI, rigidBody, quaternion, velocity, torce, cID);
}

// -------------------------------------------------------------------------------------------------
/** @brief Slaves non-master sub-body world positions and quaternions to their composite master.
    One thread per particle component (including obstacles, which are no-ops).
    @param position world-frame positions (read master, write non-masters)
    @param quaternion world-frame quaternions (read master, write non-masters)
    @param localPos per-component local position offset (body frame of composite)
    @param localQuat per-component local quaternion offset (body frame of composite)
    @param masterSlot lookup: masterSlot[compositeIdx] = current array slot of master
    @param bodyTag per-component body tag
    @param nComponents total number of components (obstacles + particles) */
template <typename T>
__GLOBAL__ void updateSubBodyPositions_Kernel(Vector3<T>*          position,
                                              Quaternion<T>*       quaternion,
                                              const Vector3<T>*    localPos,
                                              const Quaternion<T>* localQuat,
                                              const uint*          masterSlot,
                                              const uint*          bodyTag,
                                              const uint           nComponents)
{
    uint cID = blockIdx.x * blockDim.x + threadIdx.x;
    if(cID >= nComponents)
        return;

    const uint tag = bodyTag[cID];
    // Only process non-master sub-bodies
    if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
        return;

    const uint mSlot = masterSlot[getCompositeIdx(tag)];
    // world_pos  = master_pos + rotate(master_quat, local_pos)
    // world_quat = master_quat * local_quat
    position[cID]   = position[mSlot] + (quaternion[mSlot] >> localPos[cID]);
    quaternion[cID] = quaternion[mSlot] * localQuat[cID];
    T qn            = norm(quaternion[cID]);
    if(qn > EPS<T>)
        quaternion[cID] *= (T(1) / qn);
}
//@}

#endif
