#ifndef _LINKEDCELL_ATOMIC_HH_
#define _LINKEDCELL_ATOMIC_HH_

#include <cub/cub.cuh>

#include "GrainsMemBuffer.hh"
#include "LinkedCell.hh"
#include "LinkedCell_Kernels.hh"

// =================================================================================================
/** @brief The class LinkedCell_Atomic.

    This class provides functionalities to manage linked cells for collision detection in the
    simulation using an atomic approach. It is a derived class of LinkedCell and implements the
    update of linked cells based on atomic operations. This is designed to work on the device (GPU).
    This should give better time complexity than the sort-based approach.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T>
class LinkedCell_Atomic : public LinkedCell<T, MemType::DEVICE>
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
    /** \brief Temporary buffer for atomic packing: upper 32 bits = cellID, lower 32 bits =
     * particleID */
    GrainsMemBuffer<uint64_t, MemType::DEVICE> m_atomicPackBuffer;
    /** \brief Buffer of number of particles per cell */
    GrainsMemBuffer<uint, MemType::DEVICE> m_numParticlesPerCell;
    /** \brief Final cell-organized packed buffer: upper 32 bits = cellID, lower 32 bits =
     * particleID */
    GrainsMemBuffer<uint64_t, MemType::DEVICE> m_cellParticleIDs;
    /** \brief Buffer to store the prefix sums of the number of particles per cell */
    GrainsMemBuffer<uint, MemType::DEVICE> m_cellPrefixSums;
    /** \brief Buffer for atomic counters during particle writing */
    GrainsMemBuffer<uint, MemType::DEVICE> m_cellCounters;
    /** \brief Pre-allocated workspace for CUB scan operations */
    void* m_cubScanTempStorage = nullptr;
    /** \brief Size of CUB scan temporary storage */
    size_t m_cubScanTempStorageBytes = 0;
    /** \brief Device pointer for number of cells (used in resize operations) */
    uint* m_d_numCells = nullptr;
    /** \brief CUDA streams for concurrent kernel execution */
    cudaEvent_t  m_resizeComplete;  // Event to track cell resize completion
    cudaStream_t m_stream0;         // Old position copy
    cudaStream_t m_stream1;         // Neighbor cell generation
    cudaStream_t m_stream2;         // Particle packing
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    LinkedCell_Atomic() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    LinkedCell_Atomic(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rb,
                      const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                      const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& quaternions,
                      const LinkedCellParameters<T>&                         linkedCellParameters,
                      const uint                                             nObstacles,
                      const uint                                             nParticles)
        : LinkedCell<T, MemType::DEVICE>(
              rb, positions, quaternions, linkedCellParameters, nObstacles, nParticles)
        , m_atomicPackBuffer(nParticles)
        , m_numParticlesPerCell(m_numCells)
        , m_cellParticleIDs(nParticles)
        , m_cellPrefixSums(m_numCells)
        , m_cellCounters(m_numCells)
    {
        m_numParticlesPerCell.fill();

        // Pre-allocate CUB scan workspace
        // Check for errors after each batch
        cudaErrCheck(cub::DeviceScan::ExclusiveSum(nullptr,
                                                   m_cubScanTempStorageBytes,
                                                   m_numParticlesPerCell.getData(),
                                                   m_cellPrefixSums.getData(),
                                                   m_numCells));
        cudaErrCheck(cudaMalloc(&m_cubScanTempStorage, m_cubScanTempStorageBytes));

        cudaErrCheck(cudaMalloc(&m_d_numCells, sizeof(uint)));
        cudaErrCheck(cudaEventCreate(&m_resizeComplete));
        cudaErrCheck(cudaStreamCreate(&m_stream0));
        cudaErrCheck(cudaStreamCreate(&m_stream1));
        cudaErrCheck(cudaStreamCreate(&m_stream2));
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    virtual ~LinkedCell_Atomic()
    {
        if(m_cubScanTempStorage != nullptr)
            cudaErrCheck(cudaFree(m_cubScanTempStorage));
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
    /** @brief Gets number of particles per cell */
    const uint* getNumParticlesPerCell() const override
    {
        return m_numParticlesPerCell.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets packed cell-particle IDs array (upper 32 bits = cellID, lower 32 bits =
     * particleID) */
    const uint64_t* getCellParticleIDs() const
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
        +--------------------------------------+----------------------------------+
                                               \| (resize complete event)
                                               \|
        +--------------------------------------v----------------------------------+
        \| Stream 1: [wait event] -> memset(neighbors) -> generateNeighborCells      \|
        +-------------------------------------------------------------------------+

        +-------------------------------------------------------------------------+
        \| Stream 0: cudaMemcpyAsync (old positions, D2D) [independent]            \|
        +-------------------------------------------------------------------------+

        +-------------------------------------------------------------------------+
        \| Stream 2: memset(particleCounts) -> computeCellParticleIDs               \|
        +-------------------------------------------------------------------------+

        Dependencies:
        - Stream 1 waits for resize completion (cudaStreamWaitEvent)
        - Stream 2's memset must complete before pack (atomic counting needs zero)
        - Streams 0, 1, 2 run concurrently
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
            m_neighborCells.reserve(m_numCells * 27, m_stream1);
            m_numParticlesPerCell.reserve(m_numCells, m_stream2);
            m_cellPrefixSums.reserve(m_numCells, m_stream2);
            m_cellCounters.reserve(m_numCells, m_stream2);

            // Copy old positions on stream0 (independent)
            cudaMemcpyAsync(m_oldPosition.getData(),
                            m_positions->getData(),
                            (m_numObstacles + m_numParticles) * sizeof(Vector3<T>),
                            cudaMemcpyDeviceToDevice,
                            m_stream0);

            // Stream 1: Generate neighbor cells (waits for resize)
            uint numBlocksCells, numThreadsCells;
            computeOptimalThreadsAndBlocks(m_numCells,
                                           GrainsParameters<T>::m_GPU,
                                           numBlocksCells,
                                           numThreadsCells);
            cudaStreamWaitEvent(m_stream1, m_resizeComplete, 0);
            m_neighborCells.fill(UINT_MAX, m_stream1);
            generateNeighborCells_Device<<<numBlocksCells, numThreadsCells, 0, m_stream1>>>(
                m_cells.getData(),
                m_numCells,
                m_neighborCells.getData());
        }
        else
        {
            // Non-adaptive: cell geometry is fixed, no need to regenerate neighbor cells.
            // Explicitly synchronize the null stream before launching on non-null streams to
            // ensure any prior default-stream work (GJK kernels, timer syncs, etc.) has
            // completed -- this mirrors the cudaStreamSynchronize(0) that the adaptive path
            // performs inside the resize block above.
            cudaStreamSynchronize(0);
        }

        /* Launch concurrent operations */

        // Stream 2: Reset particle counts + cell counters, then pack IDs (sequential within stream)
        m_numParticlesPerCell.fill(0u, m_stream2);
        m_cellCounters.fill(0u, m_stream2);
        computeCellParticleIDs_Device<<<numBlocks, numThreads, 0, m_stream2>>>(
            m_cells.getData(),
            m_positions->getData() + m_numObstacles,
            m_numParticles,
            m_numObstacles,
            m_atomicPackBuffer.getData(),
            m_numParticlesPerCell.getData());

        // Synchronize all streams before returning
        if(m_useAdaptiveSkin)
        {
            cudaStreamSynchronize(m_stream0);  // Old position copy
            cudaStreamSynchronize(m_stream1);  // Neighbor generation
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
            // Complete preparation: resize, pack IDs, generate neighbors, link obstacles
            this->prepareLinkedCellUpdate();
            // Relink obstacles
            this->linkObstacles();
        }
        else
        {
            // Only relink obstacles if they moved
            if(SS.obstaclesMoved)
                this->linkObstacles();
            if(m_useAdaptiveSkin)
                this->launchDisplacementReduceAsync();
            return false;
        }

        // Prefix sum to find the start index of each cell in the particleIDArray
        cudaErrCheck(cub::DeviceScan::ExclusiveSum(m_cubScanTempStorage,
                                                   m_cubScanTempStorageBytes,
                                                   m_numParticlesPerCell.getData(),
                                                   m_cellPrefixSums.getData(),
                                                   m_numCells));

        // Write the particle IDs into the particleInCells using packed uint64 data
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(m_numParticles,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        writeCellParticleIDs_Kernel<<<numBlocks, numThreads>>>(m_atomicPackBuffer.getData(),
                                                               m_cellPrefixSums.getData(),
                                                               m_numParticles,
                                                               m_cellParticleIDs.getData(),
                                                               m_cellCounters.getData());

        // Launch displacement reduce async so next iteration's check costs almost nothing
        if(m_useAdaptiveSkin)
            this->launchDisplacementReduceAsync();

        return true;
    }
    //@}
};

#endif