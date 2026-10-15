#include "ComponentManager.hh"
#include "BodyTag.hh"
#include "ComponentManagerCommon.hh"
#include "ComponentManagerGPU_Kernels.hh"
#include "ForceModule.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
ComponentManager<T, M>::ComponentManager(
    GrainsMemBuffer<RigidBody<T>*, M>*              rigidBody,
    GrainsMemBuffer<uint, MemType::HOST>&&          bodyTags,
    GrainsMemBuffer<Vector3<T>, MemType::HOST>&&    position,
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>&& orientation,
    GrainsMemBuffer<Vector3<T>, MemType::HOST>&&    localPos,
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>&& localQuat,
    uint                                            nObstacles,
    uint                                            nParticles,
    uint                                            nComposites,
    uint                                            nSubBodies)
    : m_rigidBody(rigidBody)
    , m_velocity(nParticles + nObstacles)
    , m_torce(nParticles + nObstacles)
    , m_masterSlot(nComposites + 1)
    , m_counts{nObstacles, nParticles, 0u, nComposites, nSubBodies}
{
    const uint nComp = nParticles + nObstacles;
    GAssert(rigidBody->getSize() == nComp, "Rigid body size mismatch");
    const uint* hbtPtr = bodyTags.getData();  // capture before potential move

    if constexpr(M == MemType::HOST)
    {
        // pointer transfer
        m_bodyTag    = std::move(bodyTags);
        m_position   = std::move(position);
        m_quaternion = std::move(orientation);
        m_localPos   = std::move(localPos);
        m_localQuat  = std::move(localQuat);

        // Populate masterSlot (m_bodyTag is HOST-accessible)
        for(uint i = 0; i < nComp; ++i)
        {
            const uint tag = m_bodyTag[i];
            if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
                m_masterSlot[getCompositeIdx(tag)] = i;
        }

        // Derive sub-body world transforms from master positions
        updateSubBodyPositions();
    }
    else
    {
        // Upload all HOST buffers to device
        m_bodyTag.initialize(nComp);
        m_bodyTag.copyFrom(bodyTags);
        m_position.initialize(nComp);
        m_position.copyFrom(position);
        m_quaternion.initialize(nComp);
        m_quaternion.copyFrom(orientation);
        m_localPos.initialize(nComp);
        m_localPos.copyFrom(localPos);
        m_localQuat.initialize(nComp);
        m_localQuat.copyFrom(localQuat);

        // Build masterSlot on host (bodyTags still HOST-accessible; not moved), then upload
        GrainsMemBuffer<uint, MemType::HOST> hMasterSlot(nComposites + 1);
        for(uint i = 0; i < nComp; ++i)
        {
            const uint tag = bodyTags[i];
            if(isSubBody(tag) && getSubBodyLocalIdx(tag) == 0u)
                hMasterSlot[getCompositeIdx(tag)] = i;
        }
        m_masterSlot.copyFrom(hMasterSlot);
    }

    // Create CollisionDetectionModule
    m_collisionDetectionModule = std::make_unique<CollisionDetectionModule<T, M>>(
        m_rigidBody,
        m_position,
        m_quaternion,
        m_bodyTag,
        hbtPtr,
        GrainsParameters<T>::m_collisionDetection,
        m_counts.numObstacles,
        m_counts.numParticles);

    // Size pair-indexed buffers to CDM's initial pair capacity
    const size_t pairCapacity = m_collisionDetectionModule->getPairBufferSize();
    m_contactInfo.initialize(pairCapacity);
    m_pairList.initialize(pairCapacity);

    // Create ForceModule
    m_forceModule = std::make_unique<ForceModule<T, M>>(pairCapacity,
                                                        GrainsParameters<T>::m_isContactWithMemory);
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Vector3<T>, M>& ComponentManager<T, M>::getLocalPos() const
{
    return m_localPos;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Quaternion<T>, M>& ComponentManager<T, M>::getLocalQuat() const
{
    return m_localQuat;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<uint, M>& ComponentManager<T, M>::getBodyTag() const
{
    return m_bodyTag;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Vector3<T>, M>& ComponentManager<T, M>::getPosition() const
{
    return m_position;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Quaternion<T>, M>& ComponentManager<T, M>::getQuaternion() const
{
    return m_quaternion;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Kinematics<T>, M>& ComponentManager<T, M>::getVelocity() const
{
    return m_velocity;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const GrainsMemBuffer<Torce<T>, M>& ComponentManager<T, M>::getTorce() const
{
    return m_torce;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const NeighborList<T, M>* ComponentManager<T, M>::getNeighborList() const
{
    return m_collisionDetectionModule ? m_collisionDetectionModule->getNeighborList() : nullptr;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
const CollisionDetectionModule<T, M>* ComponentManager<T, M>::getCollisionDetectionModule() const
{
    return m_collisionDetectionModule.get();
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfParticles() const
{
    return m_counts.numParticles;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfObstacles() const
{
    return m_counts.numObstacles;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfComposites() const
{
    return m_counts.numComposites;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
uint ComponentManager<T, M>::getNumberOfSubBodies() const
{
    return m_counts.numSubBodies;
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::copyTo_PostProcessing(
    const std::unique_ptr<ComponentManager<T, MemType::HOST>>& other)
{
    other->setPosition(m_position);
    other->setQuaternion(m_quaternion);
    other->setVelocity(m_velocity);
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::insertParticles(const std::unique_ptr<Insertion<T>>& insertionPolicy)
{
    // Insertion is a HOST-only operation. The DEVICE branch is discarded by if constexpr so
    // that type-incompatible calls to insertionPolicy->insert() are never compiled for DEVICE.
    if constexpr(M == MemType::HOST)
    {
        insertionPolicy->insert(m_rigidBody,
                                m_position,
                                m_quaternion,
                                m_velocity,
                                GrainsParameters<T>::m_collisionDetection.linkedCellParameters,
                                m_counts.numObstacles,
                                m_counts.numParticles,
                                m_bodyTag,
                                m_localPos,
                                m_localQuat);
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::detectCollisions()
{
    m_collisionDetectionModule->run(m_rigidBody->getData(),
                                    m_position,
                                    m_quaternion,
                                    m_velocity,
                                    m_torce,
                                    m_bodyTag,
                                    m_localPos,
                                    m_localQuat,
                                    m_masterSlot,
                                    m_contactInfo,
                                    m_pairList,
                                    m_counts);
}

// -------------------------------------------------------------------------------------------------
template <typename T, MemType M>
void ComponentManager<T, M>::computeContactForces(
    const GrainsMemBuffer<ContactForceModel<T>*, M>& CF)
{
    m_forceModule->run(CF,
                       m_rigidBody,
                       m_position,
                       m_velocity,
                       m_pairList,
                       m_contactInfo,
                       m_torce,
                       m_bodyTag,
                       m_masterSlot,
                       m_counts);
}

// -------------------------------------------------------------------------------------------------
// Updates the position and quaternion of master particles (non-masters are updated by
// updateSubBodyPositions after this call).
template <typename T, MemType M>
void ComponentManager<T, M>::moveParticles(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI)
{
    if constexpr(M == MemType::HOST)
    {
        for(uint pID = m_counts.numObstacles; pID < m_counts.numObstacles + m_counts.numParticles;
            ++pID)
        {
            if(isSubBody(m_bodyTag[pID]) && getSubBodyLocalIdx(m_bodyTag[pID]) != 0u)
                continue;
            moveParticles_common(TI.getData(),
                                 m_rigidBody->getData(),
                                 m_position.getData(),
                                 m_quaternion.getData(),
                                 m_velocity.getData(),
                                 m_torce.getData(),
                                 pID);
        }
    }
    else
    {
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(m_counts.numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        moveParticles_Kernel<<<numBlocks, numThreads>>>(TI.getData(),
                                                        m_rigidBody->getData(),
                                                        m_position.getData(),
                                                        m_quaternion.getData(),
                                                        m_velocity.getData(),
                                                        m_torce.getData(),
                                                        m_bodyTag.getData(),
                                                        m_counts.numObstacles,
                                                        m_counts.numParticles);
    }
    updateSubBodyPositions();
}

// -------------------------------------------------------------------------------------------------
// Performs the second velocity half-kick (KDK Step 3; no-op for single-pass schemes).
template <typename T, MemType M>
void ComponentManager<T, M>::advanceVelocity(const GrainsMemBuffer<TimeIntegrator<T>*, M>& TI)
{
    if constexpr(M == MemType::HOST)
    {
        for(uint pID = m_counts.numObstacles; pID < m_counts.numObstacles + m_counts.numParticles;
            ++pID)
        {
            if(isSubBody(m_bodyTag[pID]) && getSubBodyLocalIdx(m_bodyTag[pID]) != 0u)
                continue;
            advanceVelocity_common(TI.getData(),
                                   m_rigidBody->getData(),
                                   m_quaternion.getData(),
                                   m_velocity.getData(),
                                   m_torce.getData(),
                                   pID);
        }
    }
    else
    {
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(m_counts.numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        advanceVelocity_Kernel<<<numBlocks, numThreads>>>(TI.getData(),
                                                          m_rigidBody->getData(),
                                                          m_quaternion.getData(),
                                                          m_velocity.getData(),
                                                          m_torce.getData(),
                                                          m_bodyTag.getData(),
                                                          m_counts.numObstacles,
                                                          m_counts.numParticles);
    }
}

// -------------------------------------------------------------------------------------------------
// Slaves non-master sub-body world transforms to their composite master.
template <typename T, MemType M>
void ComponentManager<T, M>::updateSubBodyPositions()
{
    if(m_counts.numSubBodies == 0)
        return;

    if constexpr(M == MemType::HOST)
    {
        const uint nTotal = m_counts.numObstacles + m_counts.numParticles;
        for(uint cID = m_counts.numObstacles; cID < nTotal; ++cID)
        {
            const uint tag = m_bodyTag[cID];
            if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
                continue;
            const uint mSlot  = m_masterSlot[getCompositeIdx(tag)];
            m_position[cID]   = m_position[mSlot] + (m_quaternion[mSlot] >> m_localPos[cID]);
            m_quaternion[cID] = m_quaternion[mSlot] * m_localQuat[cID];
            T qn              = norm(m_quaternion[cID]);
            if(qn > T(1e-12))
                m_quaternion[cID] *= (T(1) / qn);
        }
    }
    else
    {
        const uint nTotal = m_counts.numObstacles + m_counts.numParticles;
        uint       numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nTotal, GrainsParameters<T>::m_GPU, numBlocks, numThreads);
        updateSubBodyPositions_Kernel<<<numBlocks, numThreads>>>(m_position.getData(),
                                                                 m_quaternion.getData(),
                                                                 m_localPos.getData(),
                                                                 m_localQuat.getData(),
                                                                 m_masterSlot.getData(),
                                                                 m_bodyTag.getData(),
                                                                 nTotal);
    }
}

// --------------------------------------------------------------------------------------------------
// Explicit instantiations
template class ComponentManager<float, MemType::HOST>;
template class ComponentManager<double, MemType::HOST>;
template class ComponentManager<float, MemType::DEVICE>;
template class ComponentManager<double, MemType::DEVICE>;
