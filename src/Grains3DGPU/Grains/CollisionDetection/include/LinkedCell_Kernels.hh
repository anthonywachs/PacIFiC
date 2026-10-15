#ifndef _LINKEDCELL_KERNELS_HH_
#define _LINKEDCELL_KERNELS_HH_

#include <cooperative_groups.h>

#include "Basic.hh"
#include "Cells.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"

// =================================================================================================
/** @brief The class LinkedCell_Kernels.

    This header file contains the declarations of the various kernels used for updating the linked
    cells in the simulation.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name LinkedCell_Kernels: External Kernels */
//@{
/** @brief Resizes the cells
    @param cells pointer to the Cells object
    @param cellSize new size of the cell
    @param numCells number of cells */
template <typename T>
__GLOBAL__ void resizeCells_Device(Cells<T>** cells, const T cellSize, uint* numCells)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID > 0)
        return;

    cells[0]->resize(cellSize);
    *numCells = cells[0]->getNumCells();
}

// -------------------------------------------------------------------------------------------------
/** @brief Gets the neighbor cells array
    @param cells pointer to the Cells object
    @param numCells number of cells
    @param tr transformations */
template <typename T>
__GLOBAL__ void generateNeighborCells_Device(const Cells<T>* const* cells,
                                             const uint             numCells,
                                             uint*                  neighborCells)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numCells)
        return;
    // Generate neighbor cells for the cell with index tID
    // Each thread is responsible for one cell
    cells[0]->generateNeighborCells(neighborCells, tID, tID + 1);
}

// -------------------------------------------------------------------------------------------------
/** @brief Initialize synthetic cell prefix sums for AtomicFixed layout
    @param numCells number of cells
    @param maxParticlesPerCell maximum particles per cell
    @param cellPrefixSums output array of prefix sums (cellID * maxParticlesPerCell) */
static __GLOBAL__ void initFixedCellPrefixSums_Device(const uint numCells,
                                                      const uint maxParticlesPerCell,
                                                      uint*      cellPrefixSums)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numCells)
        return;

    cellPrefixSums[tID] = tID * maxParticlesPerCell;
}

// -------------------------------------------------------------------------------------------------
/** @brief Links obstacles to cells using support function and 1-ring expansion
    @param rb pointer to the rigid bodies (obstacles)
    @param positions world-space centers of obstacles
    @param quaternions world-space orientations of obstacles
    @param cells pointer to the Cells object (grid definition)
    @param nObstacles number of obstacles
    @param maxPerObstacle maximum number of cells an obstacle can occupy
    @param obstacleID buffer for obstacle IDs and counts (uint2)
    @param obstacleCellID buffer for cell IDs that obstacles occupy */
template <typename T>
static __GLOBAL__ void linkObstacles_Device(const RigidBody<T>* const* rb,
                                            const Vector3<T>*          positions,
                                            const Quaternion<T>*       quaternions,
                                            const Cells<T>* const*     cells,
                                            const uint                 nObstacles,
                                            const uint                 maxPerObstacle,
                                            uint2*                     obstacleID,
                                            uint*                      obstacleCellID)
{
    const uint obstacleIdx = blockIdx.x;
    if(obstacleIdx >= nObstacles)
        return;

    // Only use one thread per block for simplicity
    if(threadIdx.x != 0)
        return;

    // Lambda-like device function for support computation
    auto support = [&](const Vector3<T>& worldDirection) -> Vector3<T> {
        // Transform world direction to local coordinates using inverse rotation
        const Quaternion<T>& q              = quaternions[obstacleIdx];
        const Vector3<T>     localDirection = q << worldDirection;
        Vector3<T>           supPt          = rb[obstacleIdx]->getConvex()->support(localDirection);
        transform(q, positions[obstacleIdx], supPt);
        return supPt;
    };

    // Get grid dimensions
    const uint4 numCells = cells[0]->getNumCellsPerDirection();

    // Compute AABB by querying support in all 6 axis directions
    const Vector3<T> minExt(support(Vector3<T>(-1, 0, 0))[0],
                            support(Vector3<T>(0, -1, 0))[1],
                            support(Vector3<T>(0, 0, -1))[2]);
    const Vector3<T> maxExt(support(Vector3<T>(1, 0, 0))[0],
                            support(Vector3<T>(0, 1, 0))[1],
                            support(Vector3<T>(0, 0, 1))[2]);

    // Convert world coordinates to cell coordinates
    const uint3 minCell = cells[0]->computeCellID(minExt, false);
    const uint3 maxCell = cells[0]->computeCellID(maxExt, false);

    int minX = std::max((int)minCell.x - 1, 0);
    int maxX = std::min((int)maxCell.x + 1, (int)numCells.x - 1);
    int minY = std::max((int)minCell.y - 1, 0);
    int maxY = std::min((int)maxCell.y + 1, (int)numCells.y - 1);
    int minZ = std::max((int)minCell.z - 1, 0);
    int maxZ = std::min((int)maxCell.z + 1, (int)numCells.z - 1);

    uint cellCount = 0;
    // Offset in the obstacleCellID buffer
    const uint offset = obstacleIdx * maxPerObstacle;
    // Nested loops with 1-ring expansion; guard against cellCount >= maxPerObstacle
    // (can happen when numCells > 3375 and the cap reduces maxPerObstacle below the AABB size)
    for(int x = minX; x <= maxX && cellCount < maxPerObstacle; ++x)
    {
        for(int y = minY; y <= maxY && cellCount < maxPerObstacle; ++y)
        {
            for(int z = minZ; z <= maxZ && cellCount < maxPerObstacle; ++z)
            {
                uint cellHash = cells[0]->computeCellHash(make_uint3((uint)x, (uint)y, (uint)z));
                obstacleCellID[offset + cellCount] = cellHash;
                ++cellCount;
            }
        }
    }

    // Update the obstacle buffer
    obstacleID[obstacleIdx].x = obstacleIdx;
    obstacleID[obstacleIdx].y = cellCount;
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes the cell hash for a given point
    @param cells pointer to the Cells object
    @param positions buffer of positions
    @param numParticles number of particles
    @param cellIDs particle hash
    @param numParticlesPerCell number of particles per cell */
template <typename T>
__GLOBAL__ void computeCellID_Device(const Cells<T>* const* cells,
                                     const Vector3<T>*      positions,
                                     uint                   numParticles,
                                     uint*                  cellIDs,
                                     uint*                  numParticlesPerCell)
{
    // TODO: Load cells to shared memory if needed
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use computeDenseIndex with checkIfValid=false to avoid GAssert/__trap() on
    // out-of-grid particles. Returns UINT_MAX if position is outside the grid.
    uint c = cells[0]->computeDenseIndex(positions[tID], false);
    if(c == UINT_MAX)
        return;  // Particle is outside the grid; skip silently
    cellIDs[tID] = c;
    atomicAdd(&numParticlesPerCell[c], 1);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes packed cellID and particleID for SortBased (no particle counting)
    @param cells pointer to the Cells object
    @param positions buffer of positions
    @param numParticles number of particles
    @param particleOffset offset for particle IDs (usually numObstacles)
    @param cellParticleIDs output buffer of uint64 packed IDs */
template <typename T>
__GLOBAL__ void computeCellParticleIDs_Device(const Cells<T>* const* cells,
                                              const Vector3<T>*      positions,
                                              const uint             numParticles,
                                              const uint             particleOffset,
                                              uint64_t*              cellParticleIDs)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use computeDenseIndex with checkIfValid=false to avoid GAssert/__trap() on
    // out-of-grid particles. Returns UINT_MAX if position is outside the grid.
    uint cellID = cells[0]->computeDenseIndex(positions[tID], false);
    if(cellID == UINT_MAX)
    {
        // Mark as invalid: store UINT64_MAX so downstream kernels can detect and skip
        cellParticleIDs[tID] = UINT64_MAX;
        return;
    }
    uint particleID = tID + particleOffset;

    // Pack: upper 32 bits = cellID, lower 32 bits = particleID
    cellParticleIDs[tID] = (static_cast<uint64_t>(cellID) << 32) | particleID;
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes packed cellID and particleID for Atomic (particle counting included)
    @param cells pointer to the Cells object
    @param positions buffer of positions
    @param numParticles number of particles
    @param particleOffset offset for particle IDs (usually numObstacles)
    @param cellParticleIDs output buffer of uint64 packed IDs
    @param numParticlesPerCell atomic counter for particles per cell */
template <typename T>
__GLOBAL__ void computeCellParticleIDs_Device(const Cells<T>* const* cells,
                                              const Vector3<T>*      positions,
                                              uint                   numParticles,
                                              uint                   particleOffset,
                                              uint64_t*              cellParticleIDs,
                                              uint*                  numParticlesPerCell)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use computeDenseIndex with checkIfValid=false to avoid GAssert/__trap() on
    // out-of-grid particles. Returns UINT_MAX if position is outside the grid.
    uint cellID = cells[0]->computeDenseIndex(positions[tID], false);
    if(cellID == UINT_MAX)
    {
        cellParticleIDs[tID] = UINT64_MAX;
        return;
    }
    uint particleID = tID + particleOffset;

    // Pack data
    cellParticleIDs[tID] = (static_cast<uint64_t>(cellID) << 32) | particleID;

    // Atomic count
    atomicAdd(&numParticlesPerCell[cellID], 1);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes packed cellID and particleID for AtomicFixed (no particle counting)
    @param cells pointer to the Cells object
    @param positions buffer of positions
    @param numParticles number of particles
    @param particleOffset offset for particle IDs (usually numObstacles)
    @param maxPerCell maximum particles allowed per cell
    @param cellParticleIDs storing packed uint64 (upper 32 = cellID, lower 32 = particleID)
    @param numParticlesPerCell atomic counter for particles per cell */
template <typename T>
__GLOBAL__ void computeCellParticleIDs_Device(const Cells<T>* const* cells,
                                              const Vector3<T>*      positions,
                                              uint                   numParticles,
                                              uint                   particleOffset,
                                              uint                   maxPerCell,
                                              uint64_t*              cellParticleIDs,
                                              uint*                  numParticlesPerCell)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use computeDenseIndex with checkIfValid=false to avoid GAssert/__trap() on
    // out-of-grid particles. Returns UINT_MAX if position is outside the grid.
    uint cellID = cells[0]->computeDenseIndex(positions[tID], false);
    if(cellID == UINT_MAX)
        return;  // Particle is outside the grid; skip silently
    uint particleID = tID + particleOffset;

    // Single atomic: get slot AND increment count
    uint slot = atomicAdd(&numParticlesPerCell[cellID], 1);

    // Write packed cellID-particleID directly to 2D array if within bounds
    if(slot < maxPerCell)
    {
        uint64_t packed = (static_cast<uint64_t>(cellID) << 32) | static_cast<uint64_t>(particleID);
        cellParticleIDs[cellID * maxPerCell + slot] = packed;
    }
    // Note: Overflow silently ignored (could add error counter if needed)
}

// -------------------------------------------------------------------------------------------------
/** @brief Compute cell start indices from packed uint64 keys (upper 32 bits = cellID).
    @param cellParticleIDs Array of packed uint64 IDs
    @param numParticles Number of particles
    @param cellStart Output array to store start indices for each cell hash */
static __GLOBAL__ void computeCellStart_Kernel(const uint64_t* cellParticleIDs,
                                               const uint      numParticles,
                                               uint*           cellStart)
{
    using namespace cooperative_groups;
    thread_block           cta = this_thread_block();
    extern __shared__ uint sharedHash[];  // blockSize + 1 elements
    uint                   tid = blockIdx.x * blockDim.x + threadIdx.x;

    uint hash;
    if(tid < numParticles)
    {
        // Extract cellID from upper 32 bits
        hash                        = (uint)(cellParticleIDs[tid] >> 32);
        sharedHash[threadIdx.x + 1] = hash;
        if(tid > 0 && threadIdx.x == 0)
            sharedHash[0] = (uint)(cellParticleIDs[tid - 1] >> 32);
    }
    sync(cta);

    if(tid < numParticles)
        if(tid == 0 || hash != sharedHash[threadIdx.x])
            cellStart[hash] = tid;
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes particle IDs using packed uint64 data (for Atomic approach)
    @param atomicCellParticleIDs buffer of uint64 packed data (upper 32=cellID, lower 32=particleID)
    @param prefixSums prefix sums for each cell (starting positions)
    @param numParticles total number of particles
    @param cellParticleIDs output array where particles are written by cell
    @param cellCounters temporary counter array for atomic operations */
static __GLOBAL__ void writeCellParticleIDs_Kernel(const uint64_t* atomicCellParticleIDs,
                                                   const uint*     prefixSums,
                                                   const uint      numParticles,
                                                   uint64_t*       cellParticleIDs,
                                                   uint*           cellCounters)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;

    if(tID >= numParticles)
        return;

    // Single read gets both cellID and particleID
    uint64_t packed = atomicCellParticleIDs[tID];
    uint     cellID = static_cast<uint>(packed >> 32);

    // Skip invalid cells
    if(cellID == UINT_MAX)
        return;

    // Get the starting position for this cell from prefix sums
    uint cellStart = prefixSums[cellID];

    // Use atomic to get unique position within the cell
    uint localOffset = atomicAdd(&cellCounters[cellID], 1);

    // Write packed cellID-particleID to the computed position
    uint finalPosition             = cellStart + localOffset;
    cellParticleIDs[finalPosition] = packed;
}
//@}

#endif