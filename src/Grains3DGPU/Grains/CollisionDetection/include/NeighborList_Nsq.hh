#ifndef _NEIGHBORLIST_NSQ_HH_
#define _NEIGHBORLIST_NSQ_HH_

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "NeighborList.hh"
#include "NeighborList_Kernels.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The class NeighborList_Nsq.

    This is a derived class of NeighborList. It implements the neighbor list
    creation using an O(n^2) algorithm. This is useful for systems with a small
    number of components since we bypass LinkedCell and Bounding Volume and use
    a brute force approach.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class NeighborList_Nsq : public NeighborList<T, M>
{
    using NL = NeighborList<T, M>;
    using NL::m_pairCount;
    using NL::m_pairList;

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor */
    NeighborList_Nsq() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with number of obstacles and particles
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    NeighborList_Nsq(const uint nObstacles, const uint nParticles)
    {
        // Pre-allocate with initialPairsPerParticle hint rather than the full
        // N*(N-1)/2 budget (which can be tens of GB for large N).
        // updateNeighborList will grow the buffer to the exact size on first call.
        const auto&  LCD      = GrainsParameters<T>::m_collisionDetection.linkedCellParameters;
        const size_t allPairs = static_cast<size_t>(nObstacles) * nParticles
                                + static_cast<size_t>(nParticles) * (nParticles - 1) / 2;
        const size_t initCapacity
            = std::min(allPairs,
                       static_cast<size_t>(nParticles) * LCD.initialNumberOfPairsPerParticle
                           + static_cast<size_t>(nObstacles) * nParticles);
        m_pairList.initialize(initCapacity);
        m_pairList.fill();

        *m_pairCount = 0;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~NeighborList_Nsq() override = default;
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the neighbor list
        @param positions memory buffer of positions
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    bool updateNeighborList(GrainsMemBuffer<Vector3<T>, M>& positions,
                            const uint                      nObstacles,
                            const uint                      nParticles) final
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        // Only update at the first call
        if(SS.neighborListUpdateCount == 0)
        {
            // Grow pair buffer to exact size if the lazy initial allocation was capped.
            const size_t neededPairs = static_cast<size_t>(nObstacles) * nParticles
                                       + static_cast<size_t>(nParticles) * (nParticles - 1) / 2;
            if(m_pairList.getCapacity() < neededPairs)
                m_pairList.reserve(neededPairs);
            m_pairList.setSize(neededPairs);

            if constexpr(M == MemType::HOST || M == MemType::PINNED)
            {
                updateNeighborList_Nsq_Host(nObstacles, nParticles, m_pairList.getData());
                *m_pairCount = neededPairs;
            }
            else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
            {
                uint numBlocks, numThreads;
                computeOptimalThreadsAndBlocks(nObstacles + nParticles,
                                               GrainsParameters<T>::m_GPU,
                                               numBlocks,
                                               numThreads);
                updateNeighborList_Nsq_Device<<<numBlocks, numThreads>>>(nObstacles,
                                                                         nParticles,
                                                                         m_pairList.getData());
                *m_pairCount = neededPairs;
                cudaDeviceSynchronize();
            }
            return true;
        }
        else
            return false;
    }
    //@}
};

#endif