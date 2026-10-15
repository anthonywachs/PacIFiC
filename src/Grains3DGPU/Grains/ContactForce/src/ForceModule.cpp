#include "ForceModule.hh"
#include "BodyTag.hh"
#include "ForceModuleCommon.hh"
#include "ForceModule_Kernels.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"

// -------------------------------------------------------------------------------------------------
// Constructor
template <typename T, MemType M>
ForceModule<T, M>::ForceModule(size_t pairCapacity, bool isContactWithMemory)
{
    if constexpr(M == MemType::DEVICE)
    {
        // Compaction buffers sized to pair capacity
        m_activeIndex.initialize(pairCapacity);
        m_numActivePairsMapped.initialize(1);
        *m_numActivePairsMapped.getData() = 0u;
        const size_t cubBytes
            = queryCubSelectIfTempStorageBytes<T>(static_cast<uint>(pairCapacity),
                                                  m_activeIndex.getData(),
                                                  m_numActivePairsMapped.getDeviceData());
        m_cubSelectTempStorage.initialize(cubBytes);
    }

    if(isContactWithMemory)
    {
        uint hashCapacity = static_cast<uint>(pairCapacity / 0.7);
        m_contactTable.allocate(hashCapacity, static_cast<uint>(pairCapacity));
    }
}

// -------------------------------------------------------------------------------------------------
// Resizes GPU compaction buffers when pair buffer capacity grows; no-op on the HOST path.
template <typename T, MemType M>
void ForceModule<T, M>::resizeBuffers(size_t newPairCapacity)
{
    if constexpr(M == MemType::DEVICE)
    {
        m_activeIndex.resize(newPairCapacity);
        const size_t cubBytes
            = queryCubSelectIfTempStorageBytes<T>(static_cast<uint>(newPairCapacity),
                                                  m_activeIndex.getData(),
                                                  m_numActivePairsMapped.getDeviceData());
        m_cubSelectTempStorage.resize(cubBytes);
    }
}

// -------------------------------------------------------------------------------------------------
// Periodic mark-and-sweep
template <typename T, MemType M>
void ForceModule<T, M>::cleanupContactTable()
{
    using GP = GrainsParameters<T>;
    if(!GP::m_isContactWithMemory)
        return;

    auto& SS = GP::m_simulationState;
    if(SS.neighborListUpdateCount % 1000 == 0)
        m_contactTable.markAndSweep();
}

// -------------------------------------------------------------------------------------------------
// Accumulates non-master sub-body torces into composite masters; resets sub-body torces
template <typename T, MemType M>
void ForceModule<T, M>::assembleCompositeTorces(GrainsMemBuffer<Torce<T>, M>&         torce,
                                                const GrainsMemBuffer<Vector3<T>, M>& position,
                                                const GrainsMemBuffer<uint, M>&       bodyTag,
                                                const GrainsMemBuffer<uint, M>&       masterSlot,
                                                const ComponentCounts&                counts)
{
    if(counts.numSubBodies == 0)
        return;

    const uint nTotal = counts.numObstacles + counts.numParticles;

    if constexpr(M == MemType::HOST)
    {
        for(uint cID = counts.numObstacles; cID < nTotal; ++cID)
        {
            const uint tag = bodyTag[cID];
            if(!isSubBody(tag) || getSubBodyLocalIdx(tag) == 0u)
                continue;
            const uint       mSlot = masterSlot[getCompositeIdx(tag)];
            const Vector3<T> r     = position[cID] - position[mSlot];
            const Vector3<T> f     = torce[cID].getForce();
            const Vector3<T> tau   = torce[cID].getTorque() + (r ^ f);
            torce[mSlot].addForce(f);
            torce[mSlot].addTorque(tau);
            torce[cID].reset();
        }
    }
    else if constexpr(M == MemType::DEVICE)
    {
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(nTotal, GrainsParameters<T>::m_GPU, numBlocks, numThreads);
        assembleCompositeTorces_Kernel<<<numBlocks, numThreads>>>(torce.getData(),
                                                                  position.getData(),
                                                                  masterSlot.getData(),
                                                                  bodyTag.getData(),
                                                                  nTotal);
    }
}

// -------------------------------------------------------------------------------------------------
// Runs the complete force computation pipeline.
template <typename T, MemType M>
void ForceModule<T, M>::run(const GrainsMemBuffer<ContactForceModel<T>*, M>& CF,
                            const GrainsMemBuffer<RigidBody<T>*, M>*         rigidBody,
                            const GrainsMemBuffer<Vector3<T>, M>&            position,
                            const GrainsMemBuffer<Kinematics<T>, M>&         velocity,
                            const GrainsMemBuffer<uint2, M>&                 pairList,
                            const GrainsMemBuffer<ContactInfo<T>, M>&        contactInfo,
                            GrainsMemBuffer<Torce<T>, M>&                    torce,
                            const GrainsMemBuffer<uint, M>&                  bodyTag,
                            const GrainsMemBuffer<uint, M>&                  masterSlot,
                            const ComponentCounts&                           counts)
{
    const uint numPairs     = counts.numPairs;
    const uint numObstacles = counts.numObstacles;
    const uint numParticles = counts.numParticles;
    auto&      gt           = GrainsParameters<T>::m_fmTimer;
    gt.start(FMStage::Total);

    // 1. Periodic cleanup of contact hash table
    cleanupContactTable();

    if constexpr(M == MemType::HOST)
    {
        // 2. Compute contact forces (sequential per-pair)
        ContactMemoryView<T> contactMemory = m_contactTable.getView();
        gt.start(FMStage::ComputeForces);
        for(uint i = 0; i < numPairs; ++i)
        {
            computeContactForces_common(CF.getData(),
                                        pairList.getData(),
                                        contactInfo.getData(),
                                        position.getData(),
                                        velocity.getData(),
                                        torce.getData(),
                                        contactMemory,
                                        i);
        }
        gt.stop(FMStage::ComputeForces);

        // 3. Add external forces (gravity) to moving particles
        gt.start(FMStage::ExternalForces);
        for(uint pID = numObstacles; pID < numObstacles + numParticles; ++pID)
        {
            addExternalForces_common(GrainsParameters<T>::m_gravity,
                                     rigidBody->getData(),
                                     torce.getData(),
                                     pID);
        }
        gt.stop(FMStage::ExternalForces);
    }
    else if constexpr(M == MemType::DEVICE)
    {
        using GP = GrainsParameters<T>;

        uint numThreads, numBlocks;

        if(GP::m_useCompaction)
        {
            // 2. Compact active pair indices using CUB DeviceSelect::If
            //    The OverlapNegativeSelector predicate filters pairs in a single CUB pass;
            //    no separate flag kernel is needed.
            gt.start(FMStage::FlagAndCompact);
            const uint nActive = buildCompactActiveIndexIf(contactInfo.getData(),
                                                           numPairs,
                                                           m_activeIndex.getData(),
                                                           m_numActivePairsMapped.getDeviceData(),
                                                           m_cubSelectTempStorage.getData(),
                                                           m_cubSelectTempStorage.getSize(),
                                                           m_numActivePairsMapped.getData());
            gt.stop(FMStage::FlagAndCompact);

            if(nActive > 0)
            {
                // 4. Lazily resize per-pair intermediate buffers (indexed by original pair ID)
                m_intermediateTorceA.resize(numPairs);
                m_intermediateTorceB.resize(numPairs);

                // 5. Compute contact forces for active pairs only
                computeOptimalThreadsAndBlocks(nActive, GP::m_GPU, numBlocks, numThreads);
                ContactMemoryView<T> contactMemory = m_contactTable.getView();
                gt.start(FMStage::ComputeForces);
                computeContactForces_Kernel<<<numBlocks, numThreads>>>(
                    CF.getData(),
                    pairList.getData(),
                    contactInfo.getData(),
                    m_activeIndex.getData(),
                    position.getData(),
                    velocity.getData(),
                    m_intermediateTorceA.getData(),
                    m_intermediateTorceB.getData(),
                    contactMemory,
                    nActive);
                gt.stop(FMStage::ComputeForces);

                // 6. Reduce per-pair forces to per-particle torces using atomics
                gt.start(FMStage::ReduceTorces);
                reduceTorces_Kernel<<<numBlocks, numThreads>>>(pairList.getData(),
                                                               m_activeIndex.getData(),
                                                               m_intermediateTorceA.getData(),
                                                               m_intermediateTorceB.getData(),
                                                               torce.getData(),
                                                               nActive);
                gt.stop(FMStage::ReduceTorces);
            }
        }
        else
        {
            // No-compaction path: run force kernel over all pairs, let it check contacts internally
            if(numPairs > 0)
            {
                m_intermediateTorceA.resize(numPairs);
                m_intermediateTorceB.resize(numPairs);

                computeOptimalThreadsAndBlocks(numPairs, GP::m_GPU, numBlocks, numThreads);
                ContactMemoryView<T> contactMemory = m_contactTable.getView();

                gt.start(FMStage::ComputeForces);
                computeContactForces_AllPairs_Kernel<<<numBlocks, numThreads>>>(
                    CF.getData(),
                    pairList.getData(),
                    contactInfo.getData(),
                    position.getData(),
                    velocity.getData(),
                    m_intermediateTorceA.getData(),
                    m_intermediateTorceB.getData(),
                    contactMemory,
                    numPairs);
                gt.stop(FMStage::ComputeForces);

                gt.start(FMStage::ReduceTorces);
                reduceTorces_AllPairs_Kernel<<<numBlocks, numThreads>>>(
                    pairList.getData(),
                    m_intermediateTorceA.getData(),
                    m_intermediateTorceB.getData(),
                    torce.getData(),
                    numPairs);
                gt.stop(FMStage::ReduceTorces);
            }
        }

        // 7. Add external forces (gravity) to moving particles
        gt.start(FMStage::ExternalForces);
        computeOptimalThreadsAndBlocks(numParticles, GP::m_GPU, numBlocks, numThreads);
        addExternalForces_Kernel<<<numBlocks, numThreads>>>(GP::m_gravity[X],
                                                            GP::m_gravity[Y],
                                                            GP::m_gravity[Z],
                                                            rigidBody->getData(),
                                                            torce.getData(),
                                                            numObstacles,
                                                            numParticles);
        gt.stop(FMStage::ExternalForces);
    }

    // 8. Accumulate sub-body torces into composite masters only when needed.
    if(counts.numSubBodies > 0)
    {
        gt.start(FMStage::AssembleComposites);
        assembleCompositeTorces(torce, position, bodyTag, masterSlot, counts);
        gt.stop(FMStage::AssembleComposites);
    }
    gt.stop(FMStage::Total);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiations
template class ForceModule<float, MemType::HOST>;
template class ForceModule<double, MemType::HOST>;
template class ForceModule<float, MemType::DEVICE>;
template class ForceModule<double, MemType::DEVICE>;
