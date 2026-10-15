#ifndef _FORCEMODULECOMMON_HH_
#define _FORCEMODULECOMMON_HH_

// Uncomment to enable contact debug output (CPU path only)
// #define GRAINS_CONTACT_DEBUG

#include "Basic.hh"
#include "ContactForceModel.hh"
#include "ContactInfo.hh"
#include "ContactTable.hh"
#include "GrainsParameters.hh"
#include "Kinematics.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief ForceModule common host/device helper functions.

    Header-only inline functions shared by the HOST (CPU loop) and DEVICE (GPU kernel) paths of
    ForceModule.  Functions are templated on the scalar type T and annotated with __HOSTDEVICE__
    or __device__ as appropriate so that a single definition compiles for both paths.

    Functions moved here from ComponentManagerCommon.hh:
      - computeContactForces_common  (single-torce variant, CPU path)
      - computeContactForces_common  (dual-torce variant, GPU intermediate storage path)
      - reduceTorces_common           (__device__ only)
      - addExternalForces_common

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
/** @name ForceModule common functions between host and device */
//@{

// -------------------------------------------------------------------------------------------------
/** @brief Computes the contact forces and writes directly to per-particle torce array (CPU path).
    @param CF contact force models
    @param pairList list of pairs
    @param contactInfo contact information in the world frame
    @param position position of the components
    @param velocity kinematics of the components
    @param torce per-particle torce array (modified in-place, indexed by component ID)
    @param contactMemory view of contact memory (hash table + history data)
    @param pairID ID of the pair */
template <typename T>
__HOSTDEVICE__ static INLINE void computeContactForces_common(const ContactForceModel<T>* const* CF,
                                                              const uint2*          pairList,
                                                              const ContactInfo<T>* contactInfo,
                                                              const Vector3<T>*     position,
                                                              const Kinematics<T>*  velocity,
                                                              Torce<T>*             torce,
                                                              ContactMemoryView<T>  contactMemory,
                                                              const uint            pairID)
{
    // Contact Details
    ContactInfo<T> ci = contactInfo[pairID];
    // one load of metadata
    typename ContactInfo<T>::Snapshot snapshot = ci.getSnapshot();
    bool isContact = snapshot.overlapDistance < T(0);  // is in contact / negative distance

    // Compute the forces
    if(isContact)
    {
        const uint2 pair = pairList[pairID];
        const uint  idA  = pair.x;
        const uint  idB  = pair.y;

        // velocities of the components
        const Kinematics<T>& vA(velocity[idA]);
        const Kinematics<T>& vB(velocity[idB]);
        // geometric point of contact
        const Vector3<T>& contactPt(ci.getContactPoint());
        // relative velocity at contact point
        const Vector3<T>& relVel(vA.kinematicsAtPoint(contactPt - position[idA])
                                 - vB.kinematicsAtPoint(contactPt - position[idB]));
        // relative angular velocity
        const Vector3<T>& relAngVel(vA.getAngularComponent() - vB.getAngularComponent());

#ifdef GRAINS_CONTACT_DEBUG
        {
            static uint64_t s_contactStep = 0;
            ++s_contactStep;
            const Vector3<T>& cv   = snapshot.contactVector;
            const Vector3<T>& cp   = snapshot.contactPoint;
            const Vector3<T>& posB = position[idB];
            const Vector3<T>  omB  = vB.getAngularComponent();
            printf("[CD %lu] pair=(%u,%u) "
                   "cv=(%.4e,%.4e,%.4e) overlap=%.4e "
                   "cp=(%.4e,%.4e,%.4e) "
                   "posB=(%.4e,%.4e,%.4e) "
                   "vB=(%.4e,%.4e,%.4e) omB=(%.4e,%.4e,%.4e) "
                   "relVel=(%.4e,%.4e,%.4e)\n",
                   s_contactStep,
                   idA,
                   idB,
                   (double)cv[0],
                   (double)cv[1],
                   (double)cv[2],
                   (double)snapshot.overlapDistance,
                   (double)cp[0],
                   (double)cp[1],
                   (double)cp[2],
                   (double)posB[0],
                   (double)posB[1],
                   (double)posB[2],
                   (double)vB.getTranslationalComponent()[0],
                   (double)vB.getTranslationalComponent()[1],
                   (double)vB.getTranslationalComponent()[2],
                   (double)omB[0],
                   (double)omB[1],
                   (double)omB[2],
                   (double)relVel[0],
                   (double)relVel[1],
                   (double)relVel[2]);
        }
#endif

        // Look up or create contact history entry
        ContactHistory<T>* historyPtr = nullptr;
        if(contactMemory.m_historyData != nullptr)
        {
            uint historyIndex;
            contactMemory.findOrInsert(pair, historyIndex);
            historyPtr = &(contactMemory.m_historyData[historyIndex]);
        }

        // note that we will add torce to obstacles as well.
        uint contactForceID = snapshot.contactHash;
        CF[contactForceID]->computeForces(ci,
                                          relVel,
                                          relAngVel,
                                          position[idA],
                                          position[idB],
                                          historyPtr,
                                          torce[idA],
                                          torce[idB]);
    }
    // reset the distance so we don't compute the torce twice
    ci.setOverlapDistance(T(0));
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the contact forces (writes to intermediate per-pair storage; GPU path).
    @param CF contact force models
    @param pairList list of pairs
    @param contactInfo contact information in the world frame
    @param position position of the components
    @param velocity kinematics of the components
    @param torceA intermediate torce storage for particle A in each pair (indexed by pairID)
    @param torceB intermediate torce storage for particle B in each pair (indexed by pairID)
    @param contactMemory view of contact memory (hash table + history data)
    @param pairID ID of the pair */
template <typename T>
__HOSTDEVICE__ static INLINE void computeContactForces_common(const ContactForceModel<T>* const* CF,
                                                              const uint2*          pairList,
                                                              const ContactInfo<T>* contactInfo,
                                                              const Vector3<T>*     position,
                                                              const Kinematics<T>*  velocity,
                                                              Torce<T>*             torceA,
                                                              Torce<T>*             torceB,
                                                              ContactMemoryView<T>  contactMemory,
                                                              const uint            pairID)
{
    // Contact Details
    ContactInfo<T> ci = contactInfo[pairID];
    // one load of metadata
    typename ContactInfo<T>::Snapshot snapshot = ci.getSnapshot();
    bool isContact = snapshot.overlapDistance < T(0);  // is in contact / negative distance

    // Compute the forces
    if(isContact)  // On device path, this check is redundant.
    {
        const uint2 pair = pairList[pairID];
        const uint  idA  = pair.x;
        const uint  idB  = pair.y;

        // velocities of the components
        const Kinematics<T>& vA(velocity[idA]);
        const Kinematics<T>& vB(velocity[idB]);
        // geometric point of contact
        const Vector3<T>& contactPt(ci.getContactPoint());
        // relative velocity at contact point
        const Vector3<T>& relVel(vA.kinematicsAtPoint(contactPt - position[idA])
                                 - vB.kinematicsAtPoint(contactPt - position[idB]));
        // relative angular velocity
        const Vector3<T>& relAngVel(vA.getAngularComponent() - vB.getAngularComponent());

        // Look up or create contact history entry
        ContactHistory<T>* historyPtr = nullptr;
        if(contactMemory.m_historyData != nullptr)
        {
            uint historyIndex;
            contactMemory.findOrInsert(pair, historyIndex);
            historyPtr = &(contactMemory.m_historyData[historyIndex]);
        }

        // Reset intermediate torces before accumulating -- these per-pair slots are never
        // cleared between timesteps, and computeForces uses addForce (accumulates).
        torceA[pairID].reset();
        torceB[pairID].reset();

        // note that we will add torce to obstacles as well.
        uint contactForceID = snapshot.contactHash;
        CF[contactForceID]->computeForces(ci,
                                          relVel,
                                          relAngVel,
                                          position[idA],
                                          position[idB],
                                          historyPtr,
                                          torceA[pairID],
                                          torceB[pairID]);
    }
    // reset the distance so we don't compute the torce twice
    ci.setOverlapDistance(T(0));
}

// -------------------------------------------------------------------------------------------------
/** @brief Reduces per-pair intermediate torces to per-particle torces (DEVICE kernel helper).
    @param pairList list of pairs
    @param intermediateTorceA intermediate torce storage for particle A in each pair
    @param intermediateTorceB intermediate torce storage for particle B in each pair
    @param torce final per-particle torce array (accumulated atomically)
    @param pairID ID of the pair */
template <typename T>
__device__ static INLINE void reduceTorces_common(const uint2*    pairList,
                                                  const Torce<T>* intermediateTorceA,
                                                  const Torce<T>* intermediateTorceB,
                                                  Torce<T>*       torce,
                                                  const uint      pairID)
{
    const uint2     pair = pairList[pairID];
    const uint      idA  = pair.x;
    const uint      idB  = pair.y;
    const Torce<T>& tA   = intermediateTorceA[pairID];
    const Torce<T>& tB   = intermediateTorceB[pairID];

    // Atomically accumulate all components (6 atomics per particle)
    torce[idA].addTorceAtomic(tA);
    torce[idB].addTorceAtomic(tB);
}

// -------------------------------------------------------------------------------------------------
/** @brief Adds gravity to the component.
    @param g the gravitational acceleration vector
    @param rigidBody the rigid body of the component
    @param torce the torce acting on the component
    @param cID the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void addExternalForces_common(const Vector3<T>&          g,
                                                           const RigidBody<T>* const* rigidBody,
                                                           Torce<T>*                  torce,
                                                           const uint                 cID)
{
    const RigidBody<T>* rb   = rigidBody[cID];
    const T             mass = rb->getMass();
    // Adding the gravitational force to the torce
    torce[cID].addForce(mass * g);
}
//@}

#endif
