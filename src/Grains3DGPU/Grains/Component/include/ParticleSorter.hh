#ifndef _PARTICLESORTER_HH_
#define _PARTICLESORTER_HH_

#include <cub/cub.cuh>

#include "Basic.hh"
#include "Cells.hh"
#include "CellsFactory.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "Kinematics.hh"
#include "ParticleSorter_Kernels.hh"
#include "Quaternion.hh"
#include "Torce.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief The class ParticleSorter.

    This class sorts particles based on their Morton codes (Z-order curve) to improve memory access
    patterns and cache efficiency during collision detection. It maintains a separate Cells object
    with Morton ordering and provides methods to reorder particle data arrays accordingly.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class ParticleSorter
{
protected:
    /** @name Parameters */
    //@{
    /** \brief Cells object with Morton ordering */
    GrainsMemBuffer<Cells<T, CellOrdering::MORTON>*, M> m_cells;
    /** \brief Morton codes for each particle (Z-order curve) - CUB sort input */
    GrainsMemBuffer<uint64_t, M> m_mortonCodes;
    /** \brief Auxiliary Morton code buffer - CUB sort keys_out (discarded after sort) */
    GrainsMemBuffer<uint64_t, M> m_mortonCodesAux;
    /** \brief Index array initialised to [0..N-1] before each sort - CUB values_in */
    GrainsMemBuffer<uint, M> m_sortedIndices;
    /** \brief Sorted indices output from CUB - used by gather kernels (values_out) */
    GrainsMemBuffer<uint, M> m_sortedIndicesOut;
    /** \brief Temporary buffers for gathering sorted data */
    GrainsMemBuffer<Vector3<T>, M>    m_tempPosition;
    GrainsMemBuffer<Kinematics<T>, M> m_tempVelocity;
    GrainsMemBuffer<Quaternion<T>, M> m_tempQuaternion;
    GrainsMemBuffer<Torce<T>, M>      m_tempTorce;
    GrainsMemBuffer<uint, M>          m_tempBodyTag;
    GrainsMemBuffer<Vector3<T>, M>    m_tempLocalPos;
    GrainsMemBuffer<Quaternion<T>, M> m_tempLocalQuat;
    /** \brief CUB sort temporary storage */
    void* m_cubSortTempStorage = nullptr;
    /** \brief CUB sort temporary storage bytes */
    size_t m_cubSortTempStorageBytes = 0;
    /** \brief CUDA streams for parallel gather operations (device only) */
    cudaStream_t m_streams[7];
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    ParticleSorter() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param numObstacles total number of obstacles
        @param numParticles total number of particles */
    ParticleSorter(uint numObstacles, uint numParticles)
    {
        using GP             = GrainsParameters<T>;
        const auto& CD       = GP::m_collisionDetection;
        const auto& LCParams = CD.linkedCellParameters;
        // Create Cells object with Morton ordering
        CellsFactory<T, CellOrdering::MORTON>::template create<M>(LCParams.minCorner,
                                                                  LCParams.maxCorner,
                                                                  LCParams.cellSizeFactor,
                                                                  m_cells);

        // Allocate buffers.
        // m_sortedIndices and m_sortedIndicesOut use initialize() (not reserve()) so that
        // their logical size equals numParticles, which is required by sequence().
        m_mortonCodes.initialize(numParticles);
        m_mortonCodesAux.initialize(numParticles);
        m_sortedIndices.initialize(numParticles);
        m_sortedIndicesOut.initialize(numParticles);
        m_tempPosition.reserve(numParticles);
        m_tempVelocity.reserve(numParticles);
        m_tempQuaternion.reserve(numParticles);
        m_tempTorce.reserve(numParticles);
        m_tempBodyTag.reserve(numParticles);
        m_tempLocalPos.reserve(numParticles);
        m_tempLocalQuat.reserve(numParticles);

        if constexpr(M == MemType::DEVICE)
        {
            cudaGetLastError();
            // Initialize CUB sort workspace
            uint64_t* dummyKeys   = nullptr;
            uint*     dummyValues = nullptr;
            cudaErrCheck(cub::DeviceRadixSort::SortPairs(nullptr,
                                                         m_cubSortTempStorageBytes,
                                                         dummyKeys,
                                                         dummyKeys,
                                                         dummyValues,
                                                         dummyValues,
                                                         numParticles));
            cudaMalloc(&m_cubSortTempStorage, m_cubSortTempStorageBytes);

            // Create CUDA streams for parallel gather operations
            for(int i = 0; i < 7; ++i)
                cudaStreamCreate(&m_streams[i]);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~ParticleSorter()
    {
        // Clean up the Cells object
        if constexpr(M == MemType::HOST)
        {
            if(m_cells.getSize() > 0 && m_cells[0] != nullptr)
                delete m_cells[0];
        }
        else if constexpr(M == MemType::DEVICE)
        {
            // Free CUB workspace
            if(m_cubSortTempStorage != nullptr)
            {
                cudaFree(m_cubSortTempStorage);
                m_cubSortTempStorage = nullptr;
            }

            // Destroy CUDA streams
            for(int i = 0; i < 7; ++i)
                cudaStreamDestroy(m_streams[i]);

            // Free the Cells object on device.
            // The Cells object was created via device 'new' inside createCells_Kernel, so it
            // lives on the device heap.  Host-side cudaFree cannot free device-heap allocations
            // (it would return cudaErrorInvalidValue).  Use a device kernel instead.
            if(m_cells.getSize() > 0 && m_cells.getData() != nullptr)
            {
                deleteCells_Kernel<T, CellOrdering::MORTON><<<1, 1>>>(m_cells.getData());
                cudaDeviceSynchronize();
            }
        }
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Sorts particle data arrays based on Morton codes
        @param position particle positions (input/output)
        @param velocity particle velocities (input/output)
        @param quaternion particle orientations (input/output)
        @param torce particle forces and torques (input/output)
        @param bodyTag body tags encoding shapeId + composite membership (input/output)
        @param localPos per-component local position offsets (input/output)
        @param localQuat per-component local quaternion offsets (input/output)
        @param numObstacles number of obstacles
        @param numParticles number of particles to sort */
    void sortParticles(GrainsMemBuffer<Vector3<T>, M>&    position,
                       GrainsMemBuffer<Kinematics<T>, M>& velocity,
                       GrainsMemBuffer<Quaternion<T>, M>& quaternion,
                       GrainsMemBuffer<Torce<T>, M>&      torce,
                       GrainsMemBuffer<uint, M>&          bodyTag,
                       GrainsMemBuffer<Vector3<T>, M>&    localPos,
                       GrainsMemBuffer<Quaternion<T>, M>& localQuat,
                       uint                               numObstacles,
                       uint                               numParticles)
    {
        if constexpr(M == MemType::HOST)
        {
            // Step 1: Compute Morton codes for particles only (skip obstacles)
            uint64_t*         mortonCodes = m_mortonCodes.getData();
            const Vector3<T>* pos         = position.getData();

            for(uint i = 0; i < numParticles; ++i)
                mortonCodes[i] = m_cells[0]->computeCellHash(pos[numObstacles + i]);

            // Step 2: Initialize indices [0, 1, 2, ..., numParticles-1]
            uint* indices = m_sortedIndices.getData();
            for(uint i = 0; i < numParticles; ++i)
                indices[i] = i;

            // Step 3: Sort indices based on Morton codes
            std::sort(indices, indices + numParticles, [mortonCodes](uint a, uint b) {
                return mortonCodes[a] < mortonCodes[b];
            });

            // Step 4: Gather particle arrays according to sorted indices
            Vector3<T>*    tempPos       = m_tempPosition.getData();
            Kinematics<T>* tempVel       = m_tempVelocity.getData();
            Quaternion<T>* tempQuat      = m_tempQuaternion.getData();
            Torce<T>*      tempTorce     = m_tempTorce.getData();
            uint*          tempBodyTag   = m_tempBodyTag.getData();
            Vector3<T>*    tempLocalPos  = m_tempLocalPos.getData();
            Quaternion<T>* tempLocalQuat = m_tempLocalQuat.getData();

            const Vector3<T>*    srcPos       = position.getData() + numObstacles;
            const Kinematics<T>* srcVel       = velocity.getData() + numObstacles;
            const Quaternion<T>* srcQuat      = quaternion.getData() + numObstacles;
            const Torce<T>*      srcTorce     = torce.getData() + numObstacles;
            const uint*          srcBodyTag   = bodyTag.getData() + numObstacles;
            const Vector3<T>*    srcLocalPos  = localPos.getData() + numObstacles;
            const Quaternion<T>* srcLocalQuat = localQuat.getData() + numObstacles;

            for(uint i = 0; i < numParticles; ++i)
            {
                uint idx         = indices[i];
                tempPos[i]       = srcPos[idx];
                tempVel[i]       = srcVel[idx];
                tempQuat[i]      = srcQuat[idx];
                tempTorce[i]     = srcTorce[idx];
                tempBodyTag[i]   = srcBodyTag[idx];
                tempLocalPos[i]  = srcLocalPos[idx];
                tempLocalQuat[i] = srcLocalQuat[idx];
            }

            // Step 5: Copy sorted data back to particle arrays (skip obstacles)
            Vector3<T>*    dstPos       = position.getData() + numObstacles;
            Kinematics<T>* dstVel       = velocity.getData() + numObstacles;
            Quaternion<T>* dstQuat      = quaternion.getData() + numObstacles;
            Torce<T>*      dstTorce     = torce.getData() + numObstacles;
            uint*          dstBodyTag   = bodyTag.getData() + numObstacles;
            Vector3<T>*    dstLocalPos  = localPos.getData() + numObstacles;
            Quaternion<T>* dstLocalQuat = localQuat.getData() + numObstacles;

            for(uint i = 0; i < numParticles; ++i)
            {
                dstPos[i]       = tempPos[i];
                dstVel[i]       = tempVel[i];
                dstQuat[i]      = tempQuat[i];
                dstTorce[i]     = tempTorce[i];
                dstBodyTag[i]   = tempBodyTag[i];
                dstLocalPos[i]  = tempLocalPos[i];
                dstLocalQuat[i] = tempLocalQuat[i];
            }
        }
        else if constexpr(M == MemType::DEVICE)
        {
            using GP = GrainsParameters<T>;
            // Compute grid dimensions
            uint numBlocks, threadsPerBlock;
            computeOptimalThreadsAndBlocks(numParticles, GP::m_GPU, numBlocks, threadsPerBlock);

            // Step 1: Compute Morton codes for particles only (skip obstacles)
            computeMortonCodes_Kernel<<<numBlocks, threadsPerBlock>>>(m_cells.getData(),
                                                                      position.getData()
                                                                          + numObstacles,
                                                                      m_mortonCodes.getData(),
                                                                      numParticles);
            cudaErrCheck(cudaDeviceSynchronize());

            // Step 2: Initialize indices [0, 1, 2, ..., numParticles-1]
            m_sortedIndices.sequence();
            cudaErrCheck(cudaDeviceSynchronize());

            // Step 3: Sort indices based on Morton codes using CUB.
            // keys_in/keys_out and values_in/values_out must be non-aliasing buffers.
            // m_mortonCodesAux receives the (unused) sorted keys.
            // m_sortedIndicesOut receives the sorted values used by the gather step.
            cudaErrCheck(cub::DeviceRadixSort::SortPairs(m_cubSortTempStorage,
                                                         m_cubSortTempStorageBytes,
                                                         m_mortonCodes.getData(),
                                                         m_mortonCodesAux.getData(),
                                                         m_sortedIndices.getData(),
                                                         m_sortedIndicesOut.getData(),
                                                         numParticles));
            cudaErrCheck(cudaDeviceSynchronize());

            // Step 4: Gather particle arrays according to sorted indices using
            // multiple streams (skip obstacles).
            // m_sortedIndicesOut holds the CUB output (sorted original indices).
            // Stream 0: Position
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[0]>>>(
                position.getData() + numObstacles,
                m_tempPosition.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 1: Velocity
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[1]>>>(
                velocity.getData() + numObstacles,
                m_tempVelocity.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 2: Quaternion
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[2]>>>(
                quaternion.getData() + numObstacles,
                m_tempQuaternion.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 3: Torce
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[3]>>>(
                torce.getData() + numObstacles,
                m_tempTorce.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 4: BodyTag
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[4]>>>(
                bodyTag.getData() + numObstacles,
                m_tempBodyTag.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 5: LocalPos
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[5]>>>(
                localPos.getData() + numObstacles,
                m_tempLocalPos.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Stream 6: LocalQuat
            gather_Kernel<<<numBlocks, threadsPerBlock, 0, m_streams[6]>>>(
                localQuat.getData() + numObstacles,
                m_tempLocalQuat.getData(),
                m_sortedIndicesOut.getData(),
                numParticles);

            // Wait for all gather operations to complete
            for(int i = 0; i < 7; ++i)
                cudaStreamSynchronize(m_streams[i]);

            // Step 5: Copy sorted data back to particle arrays using multiple
            // streams (skip obstacles)
            // Stream 0: Position
            cudaMemcpyAsync(position.getData() + numObstacles,
                            m_tempPosition.getData(),
                            numParticles * sizeof(Vector3<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[0]);

            // Stream 1: Velocity
            cudaMemcpyAsync(velocity.getData() + numObstacles,
                            m_tempVelocity.getData(),
                            numParticles * sizeof(Kinematics<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[1]);

            // Stream 2: Quaternion
            cudaMemcpyAsync(quaternion.getData() + numObstacles,
                            m_tempQuaternion.getData(),
                            numParticles * sizeof(Quaternion<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[2]);

            // Stream 3: Torce
            cudaMemcpyAsync(torce.getData() + numObstacles,
                            m_tempTorce.getData(),
                            numParticles * sizeof(Torce<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[3]);

            // Stream 4: BodyTag
            cudaMemcpyAsync(bodyTag.getData() + numObstacles,
                            m_tempBodyTag.getData(),
                            numParticles * sizeof(uint),
                            cudaMemcpyDeviceToDevice,
                            m_streams[4]);

            // Stream 5: LocalPos
            cudaMemcpyAsync(localPos.getData() + numObstacles,
                            m_tempLocalPos.getData(),
                            numParticles * sizeof(Vector3<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[5]);

            // Stream 6: LocalQuat
            cudaMemcpyAsync(localQuat.getData() + numObstacles,
                            m_tempLocalQuat.getData(),
                            numParticles * sizeof(Quaternion<T>),
                            cudaMemcpyDeviceToDevice,
                            m_streams[6]);

            // Wait for all copy operations to complete
            for(int i = 0; i < 7; ++i)
                cudaStreamSynchronize(m_streams[i]);
        }
    }
    //@}
};

#endif
