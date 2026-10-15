#ifndef _RIGIDBODYFACTORY_HH_
#define _RIGIDBODYFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "RigidBody.hh"

/** @brief Unified input/output data for one category of rigid bodies (obstacles or particles).

    For obstacle and standalone-particle templates:
    - subBodiesPerTemplate = 1, isComposite = 0
    - refLocalPos/refLocalQuat = zero / identity
    - refRB, refLocalPos, refLocalQuat have size numTemplates

    For composite-particle templates:
    - subBodiesPerTemplate = S_i, isComposite = 1
    - refRB, refLocalPos, refLocalQuat are laid out flat over all sub-body prototypes;
      template i occupies indices [sum_{j<i} S_j, sum_{j<=i} S_j).

    Per-template arrays (size = numTemplates): subBodiesPerTemplate, numEachRef, isComposite,
    refInitialPos, refInitialOri.
    Per-sub-body arrays (size = sum of S_i):   refRB, refLocalPos, refLocalQuat. */
template <typename T>
struct ParticleData
{
    /** \brief Flat array of prototype rigid bodies (one per sub-body prototype). */
    GrainsMemBuffer<RigidBody<T>*> refRB;
    /** \brief Number of sub-bodies per template (1 for standalone/obstacle). */
    GrainsMemBuffer<uint> subBodiesPerTemplate;
    /** \brief Number of instances to spawn per template (Number= XML attribute). */
    GrainsMemBuffer<uint> numEachRef;
    /** \brief 1 if composite template, 0 if standalone/obstacle. */
    GrainsMemBuffer<uint> isComposite;
    /** \brief Local position offset per sub-body prototype (zero for standalone/obstacle). */
    GrainsMemBuffer<Vector3<T>> refLocalPos;
    /** \brief Local quaternion offset per sub-body prototype (identity for standalone/obstacle). */
    GrainsMemBuffer<Quaternion<T>> refLocalQuat;
    /** \brief Initial world CoM position per template. */
    GrainsMemBuffer<Vector3<T>> refInitialPos;
    /** \brief Initial world orientation per template. */
    GrainsMemBuffer<Quaternion<T>> refInitialOri;
    /** \brief Number of templates (= size of per-template arrays). */
    uint numTemplates = 0;
};

// =================================================================================================
/** @brief The class RigidBodyFactory.

    Creates the rigid body for each particle and obstacle.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T>
class RigidBodyFactory
{
private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    RigidBodyFactory() = default;

    /** @brief Destructor (forbidden) */
    ~RigidBodyFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates and returns a buffer of reference rigid bodies given an XML node
        @param obstacles XML node containing obstacle definitions
        @param particles XML node containing particle definitions
        @param obstacleData output struct to hold parsed obstacle data
        @param particleData output struct to hold parsed particle data
        @param numObstacles output number of obstacles
        @param numParticles output number of particles */
    static void create(DOMNode*         obstacles,
                       DOMNode*         particles,
                       ParticleData<T>& obstacleData,
                       ParticleData<T>& particleData,
                       uint&            numObstacles,
                       uint&            numParticles);

    /** @brief RigidBody objects must be instantiated on device, if we want to use them on device.
        Copying from host is not supported due to runtime polymorphism for this class. This
        function reads a host-side RigidBody object, and mimics it in a given device buffer. It
        calls a device kernel that is implemented in the source file.
        @param h_RB Host-side RigidBody object
        @param d_RB Device-side RigidBody object */
    static void copyHostToDevice(GrainsMemBuffer<RigidBody<T>*, MemType::HOST>&   h_RB,
                                 GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& d_RB);

    /** @brief Frees RigidBody objects that were created on device via copyHostToDevice.
        Launches a kernel that calls device-side delete on every pointer, then releases
        the pointer array.  Safe to call even if d_RB is empty.
        @param d_RB Device-side RigidBody pointer buffer */
    static void freeDevice(GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& d_RB);
    //@}
};

#endif
