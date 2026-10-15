#ifndef _NEIGHBORLIST_LINKEDCELL_HH_
#define _NEIGHBORLIST_LINKEDCELL_HH_

#include <cub/cub.cuh>

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "LinkedCell.hh"
#include "LinkedCellFactory.hh"
#include "LinkedCell_AtomicFixed.hh"
#include "LinkedCell_Host.hh"
#include "NeighborList.hh"
#include "NeighborList_Kernels.hh"

// =================================================================================================
/** @brief The class NeighborList_LinkedCell.

    This is a derived class of NeighborList. It implements the neighbor list creation using an O(n)
    algorithm. This is useful for systems with large number of components.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class NeighborList_LinkedCell : public NeighborList<T, M>
{
    using NL = NeighborList<T, M>;
    using NL::m_pairCount;
    using NL::m_pairList;

protected:
    /** @name Parameters */
    //@{
    /** \brief LinkedCell */
    std::unique_ptr<LinkedCell<T, M>> m_LinkedCell;
    /** \brief Buffer of number of neighbors for each particle */
    GrainsMemBuffer<uint, M> m_numNeighbors;
    /** \brief Buffer of prefix sums for neighbor counts */
    GrainsMemBuffer<uint, M> m_numNeighborsPrefixSums;
    /** \brief Pre-allocated workspace for CUB scan operations */
    void* m_cubScanTempStorage = nullptr;
    /** \brief Size of CUB scan temporary storage */
    size_t m_cubScanTempStorageBytes = 0;
    /** \brief Number of obstacle-particle pairs */
    uint* m_obstacleParticlePairCount = nullptr;
    /** \brief CUDA stream for obstacle-particle pair generation */
    cudaStream_t m_stream0;
    /** \brief CUDA stream for particle-particle neighbor generation */
    cudaStream_t m_stream1;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor */
    NeighborList_LinkedCell() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    NeighborList_LinkedCell(const GrainsMemBuffer<RigidBody<T>*, M>* rb,
                            const GrainsMemBuffer<Vector3<T>, M>&    positions,
                            const GrainsMemBuffer<Quaternion<T>, M>& quaternions,
                            const LinkedCellParameters<T>&           linkedCellParameters,
                            const uint                               nObstacles,
                            const uint                               nParticles)
    {
        // Create the LinkedCell buffer
        m_LinkedCell = LinkedCellFactory<T, M>::create(rb,
                                                       positions,
                                                       quaternions,
                                                       linkedCellParameters,
                                                       nObstacles,
                                                       nParticles);

        // Initialize pair list
        auto initNumPairs = GrainsParameters<T>::m_collisionDetection.linkedCellParameters
                                .initialNumberOfPairsPerParticle;
        m_pairList.initialize(initNumPairs * nParticles);
        m_pairList.fill();

        if constexpr(M == MemType::DEVICE)
        {
            // Initialize neighbor counting buffers
            m_numNeighbors.initialize(nParticles);
            m_numNeighbors.fill();

            m_numNeighborsPrefixSums.initialize(nParticles + 1);  // +1 for total
            m_numNeighborsPrefixSums.fill();
        }

        // Allocate obstacle-particle pair count
        if constexpr(M == MemType::DEVICE)
        {
            // Pre-allocate CUB scan workspace
            // First, query required workspace size (pass nullptr for temp storage)
            m_cubScanTempStorage = nullptr;
            cudaErrCheck(cub::DeviceScan::ExclusiveSum(m_cubScanTempStorage,
                                                       m_cubScanTempStorageBytes,
                                                       m_numNeighbors.getData(),
                                                       m_numNeighborsPrefixSums.getData(),
                                                       nParticles));
            // Allocate the workspace
            if(m_cubScanTempStorageBytes > 0)
            {
                cudaErrCheck(cudaMalloc(&m_cubScanTempStorage, m_cubScanTempStorageBytes));
            }

            // Allocate obstacle-particle pair count
            cudaErrCheck(cudaMallocManaged(&m_obstacleParticlePairCount, sizeof(uint)));

            // Create CUDA streams for concurrent execution
            cudaErrCheck(cudaStreamCreate(&m_stream0));
            cudaErrCheck(cudaStreamCreate(&m_stream1));
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~NeighborList_LinkedCell() override
    {
        if constexpr(M == MemType::DEVICE)
        {
            if(m_cubScanTempStorage != nullptr)
            {
                cudaErrCheck(cudaFree(m_cubScanTempStorage));
                m_cubScanTempStorage = nullptr;
            }
            if(m_obstacleParticlePairCount != nullptr)
            {
                cudaErrCheck(cudaFree(m_obstacleParticlePairCount));
                m_obstacleParticlePairCount = nullptr;
            }
            cudaErrCheck(cudaStreamDestroy(m_stream0));
            cudaErrCheck(cudaStreamDestroy(m_stream1));
        }
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the neighbor list
        @param positions array of positions
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    bool updateNeighborList(GrainsMemBuffer<Vector3<T>, M>& positions,
                            const uint                      nObstacles,
                            const uint                      nParticles) final
    {
        // Update linked cells (internally uses streams for GPU variants)
        bool LC_updated = m_LinkedCell->updateLinkedCells();

        if(LC_updated)
        {
            if constexpr(M == MemType::HOST)
            {
                auto* LC_host = static_cast<LinkedCell_Host<T>*>(m_LinkedCell.get());
                updateNeighborList_LC_Host(LC_host->getCellNeighborsList(),
                                           LC_host->getObstacleIDs(),
                                           LC_host->getObstacleCellIDs(),
                                           LC_host->getParticleIDs(),
                                           LC_host->getCellIDs(),
                                           LC_host->getCellParticles(),
                                           LC_host->getMaxCellsPerObstacle(),
                                           nObstacles,
                                           nParticles,
                                           m_pairList);
                *m_pairCount = m_pairList.getSize();
            }
            else if constexpr(M == MemType::DEVICE)
            {
                using GP = GrainsParameters<T>;
                auto& CD = GP::m_collisionDetection;
                auto& LC = CD.linkedCellParameters;

                // Reset pair count
                *m_pairCount = 0;

                // Get linked cell data
                const auto neighborCellsList = m_LinkedCell->getCellNeighborsList();
                const auto cellParticleIDs   = m_LinkedCell->getCellParticleIDs();
                const auto cellPrefixSums    = m_LinkedCell->getCellPrefixSums();
                const auto numCells          = m_LinkedCell->getNumCells();

                // Stream 0: Generate obstacle-particle pairs (concurrent with particle processing)
                if(nObstacles > 0)
                {
                    // Reset on host BEFORE kernel to avoid inter-block race: block 0's
                    // __syncthreads() is intra-block only and cannot guarantee other blocks
                    // see the zero before they call atomicAdd.
                    *m_obstacleParticlePairCount = 0;
                    // Unified obstacle-particle pair generation for all LinkedCell variants
                    generateObstacleParticlePairs_Device<<<nObstacles, 64, 0, m_stream0>>>(
                        m_LinkedCell->getObstacleIDs(),
                        m_LinkedCell->getObstacleCellIDs(),
                        cellParticleIDs,
                        cellPrefixSums,
                        m_LinkedCell->getMaxCellsPerObstacle(),
                        nObstacles,
                        nParticles,
                        numCells,
                        m_pairList.getData(),
                        m_obstacleParticlePairCount);
                    cudaStreamSynchronize(m_stream0);
                }
                else
                {
                    *m_obstacleParticlePairCount = 0;
                }

                // Stream 1: Two-phase atomic-free particle-particle neighbor generation
                uint numBlocks, numThreads;
                computeOptimalThreadsAndBlocks(nParticles,
                                               GrainsParameters<T>::m_GPU,
                                               numBlocks,
                                               numThreads);

                // Phase 1: Count neighbors per particle (on stream 1)
                m_numNeighbors.fill(0u, m_stream1);
                if(LC.type == LinkedCellType::ATOMICFIXED)
                {
                    auto* LC_atomicFixed
                        = static_cast<LinkedCell_AtomicFixed<T>*>(m_LinkedCell.get());
                    countNeighbors_AtomicFixed_Device<<<numBlocks, numThreads, 0, m_stream1>>>(
                        neighborCellsList,
                        LC_atomicFixed->getCellParticleIDs(),   // 2D fixed layout
                        LC_atomicFixed->getAtomicPackBuffer(),  // sequential layout
                        LC_atomicFixed->getNumParticlesPerCell(),
                        LC_atomicFixed->getMaxParticlesPerCell(),
                        nParticles,
                        nObstacles,
                        numCells,
                        m_numNeighbors.getData());
                }
                else
                {
                    countNeighbors_Device<<<numBlocks, numThreads, 0, m_stream1>>>(
                        neighborCellsList,
                        cellParticleIDs,
                        cellPrefixSums,
                        nParticles,
                        nObstacles,
                        numCells,
                        m_numNeighbors.getData());
                }

                // Phase 2: Compute prefix sum using CUB with pre-allocated workspace
                cudaErrCheck(cub::DeviceScan::ExclusiveSum(m_cubScanTempStorage,
                                                           m_cubScanTempStorageBytes,
                                                           m_numNeighbors.getData(),
                                                           m_numNeighborsPrefixSums.getData(),
                                                           nParticles,
                                                           m_stream1));

                // Get total pair count using async copy from exclusive scan result (on stream 1)
                // Async copy last elements: prefix_sum[n-1] + neighbor_count[n-1] = total
                uint lastPrefixSum, lastNeighborCount;
                cudaMemcpyAsync(&lastPrefixSum,
                                &m_numNeighborsPrefixSums.getData()[nParticles - 1],
                                sizeof(uint),
                                cudaMemcpyDeviceToHost,
                                m_stream1);
                cudaMemcpyAsync(&lastNeighborCount,
                                &m_numNeighbors.getData()[nParticles - 1],
                                sizeof(uint),
                                cudaMemcpyDeviceToHost,
                                m_stream1);

                // Synchronize both streams before computing total and resizing
                cudaStreamSynchronize(m_stream0);  // Wait for obstacle pairs
                cudaStreamSynchronize(m_stream1);  // Wait for prefix sum and memcpy

                uint totalPairs = lastPrefixSum + lastNeighborCount;
                // Add obstacle-particle pairs
                if(nObstacles > 0)
                    totalPairs += *m_obstacleParticlePairCount;

                // Update the pair count
                *m_pairCount = totalPairs;

                // Resize pair list to exact size (preserves existing obstacle-particle data)
                if(totalPairs != m_pairList.getSize())
                    m_pairList.resize(totalPairs);

                // Phase 3: Write neighbor pairs using prefix sums (on stream 1)
                if(LC.type == LinkedCellType::ATOMICFIXED)
                {
                    auto* LC_atomicFixed
                        = static_cast<LinkedCell_AtomicFixed<T>*>(m_LinkedCell.get());
                    updateNeighborList_LC_AtomicFixed_Device<<<numBlocks,
                                                               numThreads,
                                                               0,
                                                               m_stream1>>>(
                        neighborCellsList,
                        LC_atomicFixed->getCellParticleIDs(),   // 2D fixed layout
                        LC_atomicFixed->getAtomicPackBuffer(),  // sequential layout
                        LC_atomicFixed->getNumParticlesPerCell(),
                        m_numNeighborsPrefixSums.getData(),
                        LC_atomicFixed->getMaxParticlesPerCell(),
                        *m_obstacleParticlePairCount,
                        nObstacles,
                        nParticles,
                        numCells,
                        m_pairList.getData());
                }
                else
                {
                    updateNeighborList_LC_Device<<<numBlocks, numThreads, 0, m_stream1>>>(
                        neighborCellsList,
                        cellParticleIDs,
                        cellPrefixSums,
                        m_numNeighborsPrefixSums.getData(),
                        *m_obstacleParticlePairCount,
                        nObstacles,
                        nParticles,
                        numCells,
                        m_pairList.getData());
                }
                // Final synchronization of stream 1 (stream 0 already synchronized)
                cudaStreamSynchronize(m_stream1);
            }

            return true;
        }
        else
            return false;
    }
    //@}
};

#endif