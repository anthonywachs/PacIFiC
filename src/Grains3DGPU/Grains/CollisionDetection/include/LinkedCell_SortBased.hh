#ifndef _LINKEDCELL_SORTBASED_HH_
#define _LINKEDCELL_SORTBASED_HH_

#include <cub/cub.cuh>

#include "GrainsMemBuffer.hh"
#include "LinkedCell.hh"
#include "LinkedCell_Kernels.hh"

// =================================================================================================
/** @brief The class LinkedCell_SortBased.

    This class provides functionalities to manage linked cells for collision detection in the
    simulation using a sort-based approach. It is a derived class of LinkedCell and implements the
    update of linked cells based on sorting the particle hashes. This is designed to work on the
    device (GPU). This has optimal space complexity, while the time complexity is O(n) for computing
    the particle hashes and O(nk) for sorting the particle IDs using a parallel radix sort
    algorithm, where n is the number of particles and k is the number of cells.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T>
class LinkedCell_SortBased : public LinkedCell<T, MemType::DEVICE>
{
    using LC = LinkedCell<T, MemType::DEVICE>;
    using LC::m_cells;
    using LC::m_cellSizeWithoutSkin;
    using LC::m_maxDisplacementSquared;
    using LC::m_neighborCells;
    using LC::m_numCells;
    using LC::m_numIterationsSinceLastUpdate;
    using LC::m_numObstacles;
    using LC::m_numParticles;
    using LC::m_oldPosition;
    using LC::m_positions;
    using LC::m_skinThickness;
    using LC::m_useAdaptiveSkin;

protected:
    /** @name Parameters */
    //@{
    /** \brief Packed uint64 buffer: upper 32 bits = cellID, lower 32 bits = particleID */
    GrainsMemBuffer<uint64_t, MemType::DEVICE> m_cellParticleIDs;
    /** \brief Buffer to store cell prefix sums (start indices for each cell) */
    GrainsMemBuffer<uint, MemType::DEVICE> m_cellPrefixSums;
    /** \brief Pre-allocated workspace for CUB sort operations */
    void* m_cubSortTempStorage = nullptr;
    /** \brief Size of CUB sort temporary storage */
    size_t m_cubSortTempStorageBytes = 0;
    /** \brief Device pointer for number of cells (used in resize operations) */
    uint* m_d_numCells = nullptr;
    /** \brief CUDA streams for concurrent kernel execution */
    cudaEvent_t  m_resizeComplete;  // Event to track cell resize completion
    cudaStream_t m_stream0;         // Neighbor cell generation
    cudaStream_t m_stream1;         // Old position copy
    cudaStream_t m_stream2;         // Cell start initialization + particle packing
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    LinkedCell_SortBased() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    LinkedCell_SortBased(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rb,
                         const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                         const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& quaternions,
                         const LinkedCellParameters<T>& linkedCellParameters,
                         const uint                     nObstacles,
                         const uint                     nParticles)
        : LinkedCell<T, MemType::DEVICE>(
              rb, positions, quaternions, linkedCellParameters, nObstacles, nParticles)
        , m_cellParticleIDs(nParticles)
        , m_cellPrefixSums(m_numCells)
    {
        // Pre-allocate CUB sort workspace
        cudaErrCheck(cub::DeviceRadixSort::SortKeys(m_cubSortTempStorage,
                                                    m_cubSortTempStorageBytes,
                                                    m_cellParticleIDs.getData(),
                                                    m_cellParticleIDs.getData(),
                                                    m_numParticles));
        cudaErrCheck(cudaMalloc(&m_cubSortTempStorage, m_cubSortTempStorageBytes));

        cudaErrCheck(cudaMalloc(&m_d_numCells, sizeof(uint)));
        cudaErrCheck(cudaEventCreate(&m_resizeComplete));
        cudaErrCheck(cudaStreamCreate(&m_stream0));
        cudaErrCheck(cudaStreamCreate(&m_stream1));
        cudaErrCheck(cudaStreamCreate(&m_stream2));
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    virtual ~LinkedCell_SortBased()
    {
        if(m_cubSortTempStorage != nullptr)
            cudaErrCheck(cudaFree(m_cubSortTempStorage));
        if(m_d_numCells != nullptr)
            cudaErrCheck(cudaFree(m_d_numCells));
        cudaErrCheck(cudaEventDestroy(m_resizeComplete));
        cudaErrCheck(cudaStreamDestroy(m_stream0));
        cudaErrCheck(cudaStreamDestroy(m_stream1));
        cudaErrCheck(cudaStreamDestroy(m_stream2));
    }
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Gets packed cell-particle IDs. */
    const uint64_t* getCellParticleIDs() const override
    {
        return m_cellParticleIDs.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell prefix sums (start indices for each cell) */
    const uint* getCellPrefixSums() const override
    {
        return m_cellPrefixSums.getData();
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Prepares linked cell update: handles cell resizing, neighbor recomputation

        Concurrent execution strategy using CUDA streams:

        Timeline (adaptive skin mode):
        +-------------------------------------------------------------------------+
        \| Default Stream: resizeCells_Device -> cudaMemcpyAsync (D2H) -> [event]    \|
        +------------------------------------------+------------------------------+
                                                   \| (resize complete event)
                                                   \|
        +------------------------------------------v------------------------------+
        \| Stream 0: [wait event] -> memset(neighbors) -> generateNeighborCells      \|
        +-------------------------------------------------------------------------+

        +-------------------------------------------------------------------------+
        \| Stream 1: cudaMemcpyAsync (old positions, D2D) [independent]            \|
        +-------------------------------------------------------------------------+

        +-------------------------------------------------------------------------+
        \| Stream 2: memset(cellStart) -> computeCellParticleIDs [independent]      \|
        +-------------------------------------------------------------------------+

        Dependencies:
        - Stream 0 waits for resize completion (cudaStreamWaitEvent)
        - Streams 1 and 2 are fully independent and run concurrently
        - All streams synchronized at end before returning */
    void prepareLinkedCellUpdate()
    {
        // Compute thread configuration once
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(m_numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);

        // Full update - adjust skin thickness and resize cells
        if(m_useAdaptiveSkin)
        {
            m_skinThickness                = this->computeSkinThickness();
            T cellSize                     = m_cellSizeWithoutSkin + m_skinThickness;
            m_maxDisplacementSquared[0]    = T(0);
            m_numIterationsSinceLastUpdate = 0;

            // Resize cells on default stream (blocking but minimal impact)
            resizeCells_Device<<<1, 1>>>(m_cells.getData(), cellSize, m_d_numCells);
            cudaMemcpyAsync(&m_numCells, m_d_numCells, sizeof(uint), cudaMemcpyDeviceToHost, 0);
            cudaStreamSynchronize(0);
            cudaEventRecord(m_resizeComplete, 0);

            // Resize buffers with appropriate streams
            m_neighborCells.reserve(m_numCells * 27, m_stream0);
            m_cellPrefixSums.reserve(m_numCells, m_stream2);

            // Copy old positions on stream1 (independent)
            cudaMemcpyAsync(m_oldPosition.getData(),
                            m_positions->getData(),
                            (m_numObstacles + m_numParticles) * sizeof(Vector3<T>),
                            cudaMemcpyDeviceToDevice,
                            m_stream1);

            // Stream 0: Generate neighbor cells (waits for resize)
            uint numBlocksCells, numThreadsCells;
            computeOptimalThreadsAndBlocks(m_numCells,
                                           GrainsParameters<T>::m_GPU,
                                           numBlocksCells,
                                           numThreadsCells);
            cudaStreamWaitEvent(m_stream0, m_resizeComplete, 0);
            m_neighborCells.fill(UINT_MAX, m_stream0);
            generateNeighborCells_Device<<<numBlocksCells, numThreadsCells, 0, m_stream0>>>(
                m_cells.getData(),
                m_numCells,
                m_neighborCells.getData());
        }
        else
        {
            // Non-adaptive: cell geometry is fixed, no need to regenerate neighbor cells.
            cudaStreamSynchronize(0);
        }

        /* Launch concurrent operations */

        // Stream 2: Initialize cell start + pack particles (both independent)
        m_cellPrefixSums.fill(UINT_MAX, m_stream2);
        computeCellParticleIDs_Device<<<numBlocks, numThreads, 0, m_stream2>>>(
            m_cells.getData(),
            m_positions->getData() + m_numObstacles,
            m_numParticles,
            m_numObstacles,
            m_cellParticleIDs.getData());

        // Synchronize all streams before returning
        if(m_useAdaptiveSkin)
        {
            cudaStreamSynchronize(m_stream0);  // Neighbor generation
            cudaStreamSynchronize(m_stream1);  // Old position copy
        }
        cudaStreamSynchronize(m_stream2);  // Particle packing
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the linked cells based on the transformations */
    bool updateLinkedCells() override
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        if(this->needsUpdate())
        {
            prepareLinkedCellUpdate();
            // Relink obstacles
            this->linkObstacles();
        }
        else
        {
            // Relink obstacles if they moved
            if(SS.obstaclesMoved)
                this->linkObstacles();
            if(m_useAdaptiveSkin)
                this->launchDisplacementReduceAsync();
            return false;
        }

        // Sort packed uint64 keys using CUB with pre-allocated workspace
        cudaErrCheck(cub::DeviceRadixSort::SortKeys(m_cubSortTempStorage,
                                                    m_cubSortTempStorageBytes,
                                                    m_cellParticleIDs.getData(),
                                                    m_cellParticleIDs.getData(),
                                                    m_numParticles));
        // Compute cell start indices from sorted packed keys
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(m_numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        uint sMemSize = sizeof(uint) * (numThreads + 1);
        computeCellStart_Kernel<<<numBlocks, numThreads, sMemSize>>>(m_cellParticleIDs.getData(),
                                                                     m_numParticles,
                                                                     m_cellPrefixSums.getData());
        // Sync default stream
        cudaStreamSynchronize(0);

        // Launch displacement reduce async so next iteration's check costs almost nothing
        if(m_useAdaptiveSkin)
            this->launchDisplacementReduceAsync();

        return true;
    }
};

#endif