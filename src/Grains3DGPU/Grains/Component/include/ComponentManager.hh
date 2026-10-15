#ifndef _COMPONENTMANAGER_HH_
#define _COMPONENTMANAGER_HH_

#include "BodyTag.hh"
#include "CollisionDetectionModule.hh"
#include "ContactForceModel.hh"
#include "ContactInfo.hh"
#include "ForceModule.hh"
#include "GrainsMemBuffer.hh"
#include "Insertion.hh"
#include "Kinematics.hh"
#include "NeighborList.hh"
#include "RigidBody.hh"
#include "TimeIntegrator.hh"
#include "Torce.hh"

// =================================================================================================
/** @brief The class ComponentManager.

    Manages the state and behavior of a collection of components (particles) in the simulation.
    Each component has an associated rigid body (shape and physical properties) and per-particle
    data arrays for position, orientation, velocity, torce, etc.  The manager owns a
    CollisionDetectionModule to perform collision detection and a ForceModule to compute contact
    forces. It provides methods to initialize components, run the simulation step
    (collision detection + force computation), and update particle states.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class ComponentManager
{
protected:
    /** @name Parameters */
    //@{
    /** \brief Collision detection module */
    std::unique_ptr<CollisionDetectionModule<T, M>> m_collisionDetectionModule;
    /** \brief Force computation module */
    std::unique_ptr<ForceModule<T, M>> m_forceModule;

    /** \brief Pointer to buffer of components rigid bodies */
    const GrainsMemBuffer<RigidBody<T>*, M>* m_rigidBody;
    /** \brief Local position offset from CoM (body frame).
        For composite sub-bodies: offset from composite CoM to this sub-body CoM.
        For standalone particles: offset from user reference point to principal-axis CoM.
        Zero/identity = aligned. */
    GrainsMemBuffer<Vector3<T>, M> m_localPos;
    /** \brief Local quaternion offset (body frame).
        For composite sub-bodies: orientation of sub-body in composite body frame.
        For standalone particles: rotation from user reference frame to principal-axis frame.
        Zero/identity = aligned. */
    GrainsMemBuffer<Quaternion<T>, M> m_localQuat;
    /** \brief Body tag: encodes shape Id, composite membership and sub-body slot.
        Bit layout: [shapeId:10 | compositeIdx:14 | subBodyLocalIdx:8]
        compositeIdx==0 means standalone; use BodyTag.hh helpers to decode. */
    GrainsMemBuffer<uint, M> m_bodyTag;
    /** \brief Master slot lookup: m_masterSlot[compositeIdx] = current array slot of its master.
        Rebuilt after every Morton sort in O(N). Size = m_counts.numComposites. */
    GrainsMemBuffer<uint, M> m_masterSlot;

    /** \brief Components position */
    GrainsMemBuffer<Vector3<T>, M> m_position;
    /** \brief Components quaternion */
    GrainsMemBuffer<Quaternion<T>, M> m_quaternion;
    /** \brief Components velocities */
    GrainsMemBuffer<Kinematics<T>, M> m_velocity;
    /** \brief Components torce */
    GrainsMemBuffer<Torce<T>, M> m_torce;

    /** \brief Per-pair contact information in world frame */
    GrainsMemBuffer<ContactInfo<T>, M> m_contactInfo;
    /** \brief Pair list (indices of interacting particle pairs, populated by CDModule) */
    GrainsMemBuffer<uint2, M> m_pairList;

    /** \brief Simulation-level counts (particles, obstacles, pairs, composites, sub-bodies) */
    ComponentCounts m_counts;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    ComponentManager() = default;

    /** @brief Constructor -- fully initialises all component state and creates CDM + ForceModule.
        For HOST managers the five HOST buffers are moved (O(1) pointer swap).
        For DEVICE managers they are uploaded via cudaMemcpy.
        @param rigidBody   Pointer to the M-typed rigid body buffer
        @param bodyTags    HOST body-tag buffer (encodes shapeId); moved for HOST managers
        @param position    HOST initial position buffer; moved for HOST managers
        @param orientation HOST initial orientation buffer; moved for HOST managers
        @param localPos    HOST local-position offset buffer; moved for HOST managers
        @param localQuat   HOST local-quaternion offset buffer; moved for HOST managers
        @param nObstacles  Number of obstacles
        @param nParticles  Number of moving particles
        @param nComposites Number of composite bodies (default 0)
        @param nSubBodies  Number of sub-body slots (default 0) */
    ComponentManager(GrainsMemBuffer<RigidBody<T>*, M>*              rigidBody,
                     GrainsMemBuffer<uint, MemType::HOST>&&          bodyTags,
                     GrainsMemBuffer<Vector3<T>, MemType::HOST>&&    position,
                     GrainsMemBuffer<Quaternion<T>, MemType::HOST>&& orientation,
                     GrainsMemBuffer<Vector3<T>, MemType::HOST>&&    localPos,
                     GrainsMemBuffer<Quaternion<T>, MemType::HOST>&& localQuat,
                     uint                                            nObstacles,
                     uint                                            nParticles,
                     uint                                            nComposites = 0,
                     uint                                            nSubBodies  = 0);

    /** @brief Destructor */
    virtual ~ComponentManager() = default;
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets component local position offsets
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getLocalPos(GrainsMemBuffer<Vector3<T>, destM>& buffer) const
    {
        buffer.copyFrom(m_localPos);
    }

    /** @brief Gets component local quaternion offsets
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getLocalQuat(GrainsMemBuffer<Quaternion<T>, destM>& buffer) const
    {
        buffer.copyFrom(m_localQuat);
    }

    /** @brief Gets body tags
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getBodyTag(GrainsMemBuffer<uint, destM>& buffer) const
    {
        buffer.copyFrom(m_bodyTag);
    }

    /** @brief Gets components positions
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getPosition(GrainsMemBuffer<Vector3<T>, destM>& buffer) const
    {
        buffer.copyFrom(m_position);
    }

    /** @brief Gets components quaternions
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getQuaternion(GrainsMemBuffer<Quaternion<T>, destM>& buffer) const
    {
        buffer.copyFrom(m_quaternion);
    }

    /** @brief Gets components velocities
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getVelocity(GrainsMemBuffer<Kinematics<T>, destM>& buffer) const
    {
        m_velocity.copyTo(buffer);
    }

    /** @brief Gets components torces
        @param buffer host buffer to copy data to */
    template <MemType destM>
    void getTorce(GrainsMemBuffer<Torce<T>, destM>& buffer) const
    {
        m_torce.copyTo(buffer);
    }

    /** @brief Gets contact information in world frame
        @param buffer destination buffer to copy data into */
    template <MemType destM>
    void getContactInfo(GrainsMemBuffer<ContactInfo<T>, destM>& buffer) const
    {
        m_contactInfo.copyTo(buffer);
    }

    /** @brief Gets component local position offsets */
    const GrainsMemBuffer<Vector3<T>, M>& getLocalPos() const;

    /** @brief Gets component local quaternion offsets */
    const GrainsMemBuffer<Quaternion<T>, M>& getLocalQuat() const;

    /** @brief Gets body tags */
    const GrainsMemBuffer<uint, M>& getBodyTag() const;

    /** @brief Gets components positions */
    const GrainsMemBuffer<Vector3<T>, M>& getPosition() const;

    /** @brief Gets components quaternions */
    const GrainsMemBuffer<Quaternion<T>, M>& getQuaternion() const;

    /** @brief Gets components velocities */
    const GrainsMemBuffer<Kinematics<T>, M>& getVelocity() const;

    /** @brief Gets components torces */
    const GrainsMemBuffer<Torce<T>, M>& getTorce() const;

    /** @brief Gets neighbor list (delegated to CollisionDetectionModule) */
    const NeighborList<T, M>* getNeighborList() const;

    /** @brief Gets the collision detection module */
    const CollisionDetectionModule<T, M>* getCollisionDetectionModule() const;

    /** @brief Gets the number of particles in manager */
    uint getNumberOfParticles() const;

    /** @brief Gets the number of obstacles in manager */
    uint getNumberOfObstacles() const;

    /** @brief Gets the number of composite bodies in manager */
    uint getNumberOfComposites() const;

    /** @brief Gets the number of sub-body slots across all composites */
    uint getNumberOfSubBodies() const;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets component local position offsets
        @param p host buffer containing the local positions */
    template <MemType srcM>
    void setLocalPos(const GrainsMemBuffer<Vector3<T>, srcM>& p)
    {
        m_localPos.copyFrom(p);
    }

    /** @brief Sets component local quaternion offsets
        @param q host buffer containing the local quaternions */
    template <MemType srcM>
    void setLocalQuat(const GrainsMemBuffer<Quaternion<T>, srcM>& q)
    {
        m_localQuat.copyFrom(q);
    }

    /** @brief Sets the array of body tags
        @param id host buffer containing the body tags */
    template <MemType srcM>
    void setBodyTag(const GrainsMemBuffer<uint, srcM>& id)
    {
        m_bodyTag.copyFrom(id);
    }

    /** @brief Sets components positions
        @param p host buffer containing the positions */
    template <MemType srcM>
    void setPosition(const GrainsMemBuffer<Vector3<T>, srcM>& p)
    {
        m_position.copyFrom(p);
    }

    /** @brief Sets components quaternions
        @param q host buffer containing the quaternions */
    template <MemType srcM>
    void setQuaternion(const GrainsMemBuffer<Quaternion<T>, srcM>& q)
    {
        m_quaternion.copyFrom(q);
    }

    /** @brief Sets components velocities
        @param v host buffer containing the velocities */
    template <MemType srcM>
    void setVelocity(const GrainsMemBuffer<Kinematics<T>, srcM>& v)
    {
        m_velocity.copyFrom(v);
    }

    /** @brief Sets components torces
        @param t host buffer containing the torces */
    template <MemType srcM>
    void setTorce(const GrainsMemBuffer<Torce<T>, srcM>& t)
    {
        m_torce.copyFrom(t);
    }
    /** @brief Sets master slot lookup array
        @param ms buffer containing the master slot indices */
    template <MemType srcM>
    void setMasterSlot(const GrainsMemBuffer<uint, srcM>& ms)
    {
        m_masterSlot.copyFrom(ms);
    }
    //@}

    /** @name Manager methods */
    //@{
    /** @brief Copies particle state from this manager to another.
        @param other destination component manager */
    template <MemType srcM>
    void copyTo(const std::unique_ptr<ComponentManager<T, srcM>>& other)
    {
        other->setLocalPos(m_localPos);
        other->setLocalQuat(m_localQuat);
        other->setBodyTag(m_bodyTag);
        other->setMasterSlot(m_masterSlot);
        other->setPosition(m_position);
        other->setQuaternion(m_quaternion);
        other->setVelocity(m_velocity);
        other->setTorce(m_torce);
    }

    /** @brief Copies data to ComponentManagerCPU object for post-processing.
        @param other other component manager */
    void copyTo_PostProcessing(const std::unique_ptr<ComponentManager<T, MemType::HOST>>& other);
    //@}

    /** @name Methods */
    //@{
    /** @brief Inserts particles according to a given insertion policy
        @param ins insertion policy */
    void insertParticles(const std::unique_ptr<Insertion<T>>& insertionPolicy);

    /** @brief Runs the full collision detection pipeline via CollisionDetectionModule.
        Delegates sorting, neighbor list update, bounding volume, and narrow phase steps, as well as
        contact info transformation to world frame */
    void detectCollisions();

    /** @brief Computes contact forces and external forces via ForceModule.
        Delegates the complete force pipeline (contact forces + gravity) to m_forceModule->run().
        @param CF array of all contact force models */
    void computeContactForces(const GrainsMemBuffer<ContactForceModel<T>*, M>& CF);

    /** @brief Updates the position and velocities of particles
        @param TI time integration scheme */
    void moveParticles(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI);

    /** @brief Performs the second velocity half-kick for split-step schemes (e.g. Leapfrog).
        For single-pass schemes (e.g. FirstOrderExplicit) this is a no-op because
        TimeIntegrator::AdvanceVelocity defaults to an empty body.
        @param TI time integration scheme */
    void advanceVelocity(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI);

    /** @brief Slaves non-master sub-body positions/quaternions to their composite master.
        Must be called after moveParticles() and after advanceVelocity(). No-op if no composites. */
    void updateSubBodyPositions();
    //@}
};

#endif