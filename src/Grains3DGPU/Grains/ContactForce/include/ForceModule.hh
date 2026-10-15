#ifndef _FORCEMODULE_HH_
#define _FORCEMODULE_HH_

#include "BodyTag.hh"
#include "ContactForceModel.hh"
#include "ContactInfo.hh"
#include "ContactTable.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "Kinematics.hh"
#include "RigidBody.hh"
#include "Torce.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief The class ForceModule.

    Encapsulates the full force computation pipeline:
      1. Periodic contact-table cleanup (mark-and-sweep)
      2. [GPU only] Build compact list of active pair indices
      3. Contact force computation (parallel per-pair on GPU, sequential loop on CPU)
      4. [GPU only] Atomic reduction of per-pair intermediate torces to per-particle torces
      5. External forces (gravity applied to each moving particle)

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class ForceModule
{
private:
    /** @name Owned resources */
    //@{
    /** \brief Compact array of active pair indices (DEVICE path only) */
    GrainsMemBuffer<uint, M> m_activeIndex;
    /** \brief CUB DeviceSelect::If temporary storage (DEVICE path only) */
    GrainsMemBuffer<uint8_t, M> m_cubSelectTempStorage;
    /** \brief Zero-copy mapped pinned buffer holding the CUB compaction count.
        Device writes via getDeviceData() (mapped alias); host reads via getData(). */
    GrainsMemBuffer<uint, MemType::MAPPED> m_numActivePairsMapped;
    /** \brief Per-pair intermediate torce for particle A (lazy resize; DEVICE path only) */
    GrainsMemBuffer<Torce<T>, M> m_intermediateTorceA;
    /** \brief Per-pair intermediate torce for particle B (lazy resize; DEVICE path only) */
    GrainsMemBuffer<Torce<T>, M> m_intermediateTorceB;
    /** \brief Contact history hash table (allocated only when isContactWithMemory == true) */
    ContactHashTable<T, M> m_contactTable;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    ForceModule() = delete;

    /** @brief Constructor
        @param pairCapacity      Initial pair buffer capacity (from CDModule::getPairBufferSize())
        @param isContactWithMemory Whether contact history tracking is needed */
    ForceModule(size_t pairCapacity, bool isContactWithMemory);

    /** @brief Destructor */
    ~ForceModule() = default;

    /** @brief Deleted copy constructor */
    ForceModule(const ForceModule&) = delete;

    /** @brief Deleted copy assignment operator */
    ForceModule& operator=(const ForceModule&) = delete;

    /** @brief Defaulted move constructor */
    ForceModule(ForceModule&&) = default;

    /** @brief Defaulted move assignment operator */
    ForceModule& operator=(ForceModule&&) = default;
    //@}

    /** @name Methods */
    //@{
    /** @brief Performs periodic mark-and-sweep cleanup of the contact hash table.
        No-op when isContactWithMemory == false or the cleanup interval has not elapsed. */
    void cleanupContactTable();

    /** @brief Resizes internal GPU compaction buffers when the pair buffer capacity grows.
        @param newPairCapacity New pair buffer capacity */
    void resizeBuffers(size_t newPairCapacity);

    /** @brief Accumulates forces/torques from non-master sub-bodies into their composite master
        and resets the sub-body torces. No-op when counts.numSubBodies == 0.
        @param torce        Per-component torce array (modified in-place)
        @param position     Per-component position array
        @param bodyTag      Per-component body tag
        @param masterSlot   Per-composite master slot lookup
        @param counts       Component counts */
    void assembleCompositeTorces(GrainsMemBuffer<Torce<T>, M>&         torce,
                                 const GrainsMemBuffer<Vector3<T>, M>& position,
                                 const GrainsMemBuffer<uint, M>&       bodyTag,
                                 const GrainsMemBuffer<uint, M>&       masterSlot,
                                 const ComponentCounts&                counts);

    /** @brief Runs the complete force computation pipeline. Steps:
          1. cleanupContactTable (periodic mark-and-sweep every 1000 NL updates)
          2. [DEVICE] resize m_intermediateTorceA/B to numPairs
          3. [DEVICE] computeContactForces_Kernel -> reduceTorces_Kernel
             [HOST]   sequential computeContactForces_common loop
          4. addExternalForces (gravity)
          5. assembleCompositeTorces (no-op when numSubBodies == 0)
        @param CF            Array of contact force models
        @param rigidBody     Per-component rigid body pointer array
        @param position      Per-component position array
        @param velocity      Per-component kinematics array
        @param pairList      Per-pair component index pairs (from CDModule)
        @param contactInfo   Per-pair contact information in world frame (from CDModule)
        @param torce         Per-component torce array (modified in-place)
        @param bodyTag       Per-component body tag (encodes composite membership)
        @param masterSlot    Per-composite master slot lookup (size = numComposites)
        @param counts        Component counts (numPairs, numSubBodies, numObstacles, numParticles)
     */
    void run(const GrainsMemBuffer<ContactForceModel<T>*, M>& CF,
             const GrainsMemBuffer<RigidBody<T>*, M>*         rigidBody,
             const GrainsMemBuffer<Vector3<T>, M>&            position,
             const GrainsMemBuffer<Kinematics<T>, M>&         velocity,
             const GrainsMemBuffer<uint2, M>&                 pairList,
             const GrainsMemBuffer<ContactInfo<T>, M>&        contactInfo,
             GrainsMemBuffer<Torce<T>, M>&                    torce,
             const GrainsMemBuffer<uint, M>&                  bodyTag,
             const GrainsMemBuffer<uint, M>&                  masterSlot,
             const ComponentCounts&                           counts);
    //@}
};

#endif
