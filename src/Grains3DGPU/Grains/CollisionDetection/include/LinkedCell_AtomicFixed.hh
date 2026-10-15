#ifndef _LINKEDCELL_ATOMICFIXED_HH_
#define _LINKEDCELL_ATOMICFIXED_HH_

#include "GrainsMemBuffer.hh"
#include "LinkedCell.hh"
#include "LinkedCell_Kernels.hh"

// =================================================================================================
/** @brief The class LinkedCell_AtomicFixed.

    This class provides functionalities to manage linked cells for collision detection using a
    pre-sized fixed array approach. Unlike LinkedCell_Atomic which uses dynamic sizing with prefix
    sums, this variant allocates a fixed 2D array [numCells x maxParticlesPerCell] for direct
    particle insertion.

    Best for: Dense/uniform particle distributions where memory overhead is acceptable in exchange
    for faster updates (no prefix sum needed).

    Trade-off: Wastes memory for sparse distributions but provides O(1) insertion without prefix sum
    overhead.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T>
class LinkedCell_AtomicFixed : public LinkedCell<T, MemType::DEVICE>
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
    /** \brief Buffer of number of particles per cell */
    GrainsMemBuffer<uint, MemType::DEVICE> m_numParticlesPerCell;
    /** \brief Buffer storing packed cellID-particleID pairs: [numCells x maxParticlesPerCell]
     * (uint64: upper 32 = cellID, lower 32 = particleID) */
    GrainsMemBuffer<uint64_t, MemType::DEVICE> m_cellParticleIDs;
    /** \brief Synthetic prefix sums for uniform interface (cellID * maxParticlesPerCell) */
    GrainsMemBuffer<uint, MemType::DEVICE> m_cellPrefixSums;
    /** \brief Sequential packed buffer: atomicPackBuffer[tID] = (cellID<<32 | particleID)
     * Used by NL kernels to map thread tID -> current particle (the 2D cellParticleIDs
     * array cannot be indexed by tID since its layout is [cellID*maxPerCell+slot]). */
    GrainsMemBuffer<uint64_t, MemType::DEVICE> m_atomicPackBuffer;
    /** \brief Maximum number of particles that can fit in a single cell */
    uint m_maxParticlesPerCell;
    /** \brief Device pointer for number of cells (used in resize operations) */
    uint* m_d_numCells = nullptr;
    /** \brief CUDA streams for concurrent kernel execution */
    cudaEvent_t  m_resizeComplete;  // Event to track cell resize completion
    cudaStream_t m_stream0;         // Old position copy
    cudaStream_t m_stream1;         // Neighbor cell generation
    cudaStream_t m_stream2;         // Particle insertion
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    LinkedCell_AtomicFixed() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    LinkedCell_AtomicFixed(const GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>* rb,
                           const GrainsMemBuffer<Vector3<T>, MemType::DEVICE>&    positions,
                           const GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>& quaternions,
                           const LinkedCellParameters<T>& linkedCellParameters,
                           const uint                     nObstacles,
                           const uint                     nParticles)
        : LinkedCell<T, MemType::DEVICE>(
              rb, positions, quaternions, linkedCellParameters, nObstacles, nParticles)
        , m_numParticlesPerCell(m_numCells)
        , m_cellParticleIDs(m_numCells * linkedCellParameters.maxParticlesPerCell)
        , m_cellPrefixSums(m_numCells)
        , m_atomicPackBuffer(nParticles)
        , m_maxParticlesPerCell(linkedCellParameters.maxParticlesPerCell)
    {
        m_numParticlesPerCell.fill();
        m_cellParticleIDs.fill(UINT64_MAX);

        // Initialize synthetic prefix sums once for maximum cells (deterministic formula)
        // Since we start with smallest cell size = maximum number of cells, this covers all cases
        uint numBlocks, numThreads;
        computeOptimalThreadsAndBlocks(m_numCells,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);
        initFixedCellPrefixSums_Device<<<numBlocks, numThreads>>>(m_numCells,
                                                                  m_maxParticlesPerCell,
                                                                  m_cellPrefixSums.getData());
        cudaDeviceSynchronize();

        cudaErrCheck(cudaMalloc(&m_d_numCells, sizeof(uint)));
        cudaErrCheck(cudaEventCreate(&m_resizeComplete));
        cudaErrCheck(cudaStreamCreate(&m_stream0));
        cudaErrCheck(cudaStreamCreate(&m_stream1));
        cudaErrCheck(cudaStreamCreate(&m_stream2));
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    virtual ~LinkedCell_AtomicFixed()
    {
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
    /** @brief Gets packed cell-particle IDs array (upper 32 bits = cellID, lower 32 bits =
     * particleID) */
    const uint64_t* getCellParticleIDs() const
    {
        return m_cellParticleIDs.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets number of particles per cell */
    const uint* getNumParticlesPerCell() const override
    {
        return m_numParticlesPerCell.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets maximum particles per cell */
    uint getMaxParticlesPerCell() const
    {
        return m_maxParticlesPerCell;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets synthetic cell prefix sums (cellID * maxParticlesPerCell for fixed layout) */
    const uint* getCellPrefixSums() const override
    {
        return m_cellPrefixSums.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets sequential packed buffer (atomicPackBuffer[tID] = cellID<<32 | particleID) */
    const uint64_t* getAtomicPackBuffer() const
    {
        return m_atomicPackBuffer.getData();
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
        \| Stream 2: memset(particleCounts) -> memset(cellParticles) -> directInsert \|
        +-------------------------------------------------------------------------+

        Dependencies:
        - Stream 1 waits for resize completion (cudaStreamWaitEvent)
        - Stream 2's memsets must complete before insert (atomic insertion needs clean state)
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
            m_maxDisplacementSquared[0]    = T(0);
            m_numIterationsSinceLastUpdate = 0;

            T cellSize = m_cellSizeWithoutSkin + m_skinThickness;

            // Resize cells on default stream (blocking but minimal impact)
            resizeCells_Device<<<1, 1>>>(m_cells.getData(), cellSize, m_d_numCells);
            cudaMemcpyAsync(&m_numCells, m_d_numCells, sizeof(uint), cudaMemcpyDeviceToHost, 0);
            cudaStreamSynchronize(0);
            cudaEventRecord(m_resizeComplete, 0);

            // Resize buffers with appropriate streams
            m_neighborCells.reserve(m_numCells * 27, m_stream1);
            m_numParticlesPerCell.reserve(m_numCells, m_stream2);
            m_cellParticleIDs.reserve(m_numCells * m_maxParticlesPerCell, m_stream2);

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
            cudaStreamSynchronize(0);
        }

        // Stream 2: Reset counts + clear array + insert particles (sequential within stream)
        m_numParticlesPerCell.fill(0u, m_stream2);
        m_cellParticleIDs.fill(UINT64_MAX, m_stream2);
        computeCellParticleIDs_Device<<<numBlocks, numThreads, 0, m_stream2>>>(
            m_cells.getData(),
            m_positions->getData() + m_numObstacles,
            m_numParticles,
            m_numObstacles,
            m_maxParticlesPerCell,
            m_cellParticleIDs.getData(),
            m_numParticlesPerCell.getData());
        // Sequential layout: atomicPackBuffer[tID] = (cellID << 32 | particleID)
        // Uses the SortBased variant (no atomic counting) so tID == particle index.
        computeCellParticleIDs_Device<<<numBlocks, numThreads, 0, m_stream2>>>(
            m_cells.getData(),
            m_positions->getData() + m_numObstacles,
            m_numParticles,
            m_numObstacles,
            m_atomicPackBuffer.getData());

        // Synchronize all streams before returning
        if(m_useAdaptiveSkin)
        {
            cudaStreamSynchronize(m_stream0);  // Old position copy
            cudaStreamSynchronize(m_stream1);  // Neighbor generation
        }
        cudaStreamSynchronize(m_stream2);  // Particle insertion
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the linked cells based on the transformations */
    bool updateLinkedCells() override
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        if(this->needsUpdate())
        {
            // Complete preparation: resize, insert particles, generate neighbors, link obstacles
            this->prepareLinkedCellUpdate();
            this->linkObstacles();
            if(m_useAdaptiveSkin)
                this->launchDisplacementReduceAsync();
            return true;
        }

        // No LC rebuild needed: only relink obstacles if they moved
        if(SS.obstaclesMoved)
            this->linkObstacles();
        if(m_useAdaptiveSkin)
            this->launchDisplacementReduceAsync();
        return false;
    }
    //@}
};

#endif
