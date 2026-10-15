#include <algorithm>
#include <set>

#include "NeighborList_Kernels.hh"
#include "Transform3.hh"

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on host using an O(n^2) algorithm
__HOST__ void
    updateNeighborList_Nsq_Host(const uint nObstacles, const uint nParticles, uint2* pairList)
{
    for(uint i = 0; i < nObstacles; ++i)
        for(uint j = 0; j < nParticles; ++j)
            pairList[nParticles * i + j] = make_uint2(i, nObstacles + j);

    // Offset for p-p interactions.
    uint offset = nObstacles * nParticles;
    for(uint i = 0; i < nParticles; ++i)
        for(uint j = i + 1; j < nParticles; ++j)
            pairList[offset + i + j * (j - 1) / 2] = make_uint2(nObstacles + i, nObstacles + j);
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on device using an O(n^2) algorithm
__GLOBAL__ void
    updateNeighborList_Nsq_Device(const uint nObstacles, const uint nParticles, uint2* pairList)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID < nObstacles)
    {
        // Obstacle to particle pairs
        for(uint j = 0; j < nParticles; ++j)
            pairList[nParticles * tID + j] = make_uint2(tID, nObstacles + j);
    }
    else if(tID < nObstacles + nParticles)
    {
        // offset
        const uint offset = nObstacles * nParticles;
        // adjust tID to start from 0 for particles
        tID -= nObstacles;
        // Particle to obstacle pairs
        for(uint j = tID + 1; j < nParticles; ++j)
            pairList[offset + tID + j * (j - 1) / 2] = make_uint2(nObstacles + tID, nObstacles + j);
    }
    else
        return;
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on host using a linked cell approach
__HOST__ void updateNeighborList_LC_Host(const uint*                            cellNeighborsList,
                                         const uint2*                           obstacleIDs,
                                         const uint*                            obstacleCellIDs,
                                         const uint*                            particleIDs,
                                         const uint*                            cellIDs,
                                         const std::vector<std::list<uint>>&    cellParticles,
                                         const uint                             maxCellsPerObstacle,
                                         const uint                             numObstacles,
                                         const uint                             numParticles,
                                         GrainsMemBuffer<uint2, MemType::HOST>& pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;  // Number of neighboring cells

    // Set pairList size to zero
    pairList.setSize(0);

    // FIRST PASS: Loop over all obstacles
    for(uint i = 0; i < numObstacles; ++i)
    {
        const uint offset = i * maxCellsPerObstacle;

        const uint obstacleIndex      = obstacleIDs[i].x;
        const uint numCellsToTraverse = obstacleIDs[i].y;

        for(uint c = 0; c < numCellsToTraverse; ++c)
        {
            const uint  cell                = obstacleCellIDs[offset + c];
            const auto& targetCellParticles = cellParticles[cell];

            // Check against all particles in the target cell
            for(uint particleID : targetCellParticles)
                pairList.push_back(make_uint2(obstacleIndex, particleID));
        }
    }

    // SECOND PASS: Cell-centric approach for particle-particle pairs
    for(uint c = 0; c < cellParticles.size(); ++c)
    {
        if(cellParticles[c].empty())
            continue;  // Skip empty cells

        // Get neighbor cells for this cell (includes own cell)
        const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * c];

        // Check interactions with neighboring cells
        for(uint cc = 0; cc < NUM_NEIGHBOR_CELLS; ++cc)
        {
            uint targetCell = neighborCells[cc];
            if(targetCell == UINT_MAX || targetCell < c)
                continue;  // Skip invalid cells and cells with lower indices

            const auto& neighborCellParticles = cellParticles[targetCell];
            if(neighborCellParticles.empty())
                continue;  // Skip empty target cells

            // Process particle pairs between cells
            for(uint primaryParticle : cellParticles[c])
            {
                for(uint otherParticle : neighborCellParticles)
                {
                    if(targetCell == c)
                    {
                        // Same-cell: avoid duplicates with ordering check
                        if(primaryParticle < otherParticle)
                            pairList.push_back(make_uint2(primaryParticle, otherParticle));
                    }
                    else
                    {
                        // Cross-cell: always emit -- uniqueness is guaranteed by targetCell > c
                        const uint lo
                            = primaryParticle < otherParticle ? primaryParticle : otherParticle;
                        const uint hi
                            = primaryParticle < otherParticle ? otherParticle : primaryParticle;
                        pairList.push_back(make_uint2(lo, hi));
                    }
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Generate obstacle-particle pairs on device
__GLOBAL__ void generateObstacleParticlePairs_Device(const uint2*    obstacleIDs,
                                                     const uint*     obstacleCellIDs,
                                                     const uint64_t* packedCellParticleIDs,
                                                     const uint*     cellPrefixSums,
                                                     const uint      maxCellsPerObstacle,
                                                     const uint      numObstacles,
                                                     const uint      numParticles,
                                                     const uint      numCells,
                                                     uint2*          pairList,
                                                     uint*           pairCount)
{
    uint obstacleIdx = blockIdx.x;
    if(obstacleIdx >= numObstacles)
        return;

    const uint offset             = obstacleIdx * maxCellsPerObstacle;
    const uint obstacleIndex      = obstacleIDs[obstacleIdx].x;
    const uint numCellsToTraverse = obstacleIDs[obstacleIdx].y;

    // Each thread handles one cell for this obstacle
    for(uint c = threadIdx.x; c < numCellsToTraverse; c += blockDim.x)
    {
        const uint cell      = obstacleCellIDs[offset + c];
        const uint cellStart = cellPrefixSums[cell];

        if(cellStart == UINT_MAX)
            continue;  // Empty cell

        // Find cell end (handles both variants: with/without UINT_MAX sentinels)
        uint cellEnd;
        uint k = cell;
        do
        {
            ++k;
            cellEnd = (k < numCells) ? cellPrefixSums[k] : numParticles;
        } while(cellEnd == UINT_MAX && k < numCells);

        // Add pairs for all particles in this cell
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            // Unpack particleID from lower 32 bits
            uint64_t packed = packedCellParticleIDs[p];
            // Skip empty slots: AtomicFixed layout fills unused entries with UINT64_MAX.
            // Without this check, all maxParticlesPerCell slots are written as pairs
            // (most with particleID = UINT_MAX), causing a massive OOB overwrite of the
            // pair list buffer.
            if(packed == UINT64_MAX)
                continue;
            uint particleID       = static_cast<uint>(packed & 0xFFFFFFFF);
            uint globalIndex      = atomicAdd(pairCount, 1);
            pairList[globalIndex] = make_uint2(obstacleIndex, particleID);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Count neighbors per particle using linked cells with packed keys
__GLOBAL__ void countNeighbors_Device(const uint*     cellNeighborsList,
                                      const uint64_t* cellParticleIDs,
                                      const uint*     cellPrefixSums,
                                      const uint      numParticles,
                                      const uint      numObstacles,
                                      const uint      numCells,
                                      uint*           neighborCounts)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Unpack cellID and particleID from uint64
    const uint64_t packed = cellParticleIDs[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);  // global particle ID

    if(cell == UINT_MAX)
        return;  // neighborCounts pre-filled with 0; no write needed

    const uint* neighborCells  = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        totalNeighbors = 0;

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Calculate number of particles in cell from cellPrefixSums
        uint cellStart = cellPrefixSums[c];
        if(cellStart == UINT_MAX)
            continue;

        // Find end of cell
        uint k = c;
        uint cellEnd;
        do
        {
            ++k;
            cellEnd = cellPrefixSums[k];
        } while(cellEnd == UINT_MAX && k < numCells);

        if(k == numCells)
            cellEnd = numParticles;

        // Count particles in this neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint j = (uint)(cellParticleIDs[p] & 0xFFFFFFFF);

            if(c == cell && j <= i)
                continue;

            totalNeighbors++;
        }
    }

    neighborCounts[i - numObstacles] = totalNeighbors;
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list on device
__GLOBAL__ void updateNeighborList_LC_Device(const uint*     cellNeighborsList,
                                             const uint64_t* cellParticleIDs,
                                             const uint*     cellPrefixSums,
                                             const uint*     numNeighborsPrefixSums,
                                             const uint      pairListOffset,
                                             const uint      numObstacles,
                                             const uint      numParticles,
                                             const uint      numCells,
                                             uint2*          pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Unpack cellID and particleID from uint64
    const uint64_t packed = cellParticleIDs[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        insertIndex   = pairListOffset + numNeighborsPrefixSums[i - numObstacles];

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        uint cellStart = cellPrefixSums[c];
        if(cellStart == UINT_MAX)
            continue;

        // Find end of cell
        uint k = c;
        uint cellEnd;
        do
        {
            ++k;
            cellEnd = cellPrefixSums[k];
        } while(cellEnd == UINT_MAX && k < numCells);

        if(k == numCells)
            cellEnd = numParticles;

        // Loop through particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            // Unpack particle ID from cellParticleIDs
            const uint j = (uint)(cellParticleIDs[p] & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            // Always write (min,max) to keep pairs canonically ordered.
            pairList[insertIndex++] = (i < j) ? make_uint2(i, j) : make_uint2(j, i);
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Count neighbors per particle for AtomicFixed using direct count lookup
__GLOBAL__ void countNeighbors_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                  const uint64_t* cellParticleIDsFixed,
                                                  const uint64_t* particleCellIDsSeq,
                                                  const uint*     numParticlesPerCell,
                                                  const uint      maxParticlesPerCell,
                                                  const uint      numParticles,
                                                  const uint      numObstacles,
                                                  const uint      numCells,
                                                  uint*           neighborCounts)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use sequential buffer: particleCellIDsSeq[tID] = (cellID<<32 | particleID)
    // cellParticleIDsFixed uses 2D layout [cellID*maxPerCell+slot] for neighbour traversal.
    const uint64_t packed = particleCellIDsSeq[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        count         = 0;

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Compute cellStart for fixed layout: cellID * maxParticlesPerCell
        uint cellStart = c * maxParticlesPerCell;
        uint cellEnd   = cellStart + numParticlesPerCell[c];

        // Count particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint j = (uint)(cellParticleIDsFixed[p] & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            count++;
        }
    }

    neighborCounts[i - numObstacles] = count;
}

// -------------------------------------------------------------------------------------------------
// Updates the neighbor list for AtomicFixed using direct count lookup
__GLOBAL__ void updateNeighborList_LC_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                         const uint64_t* cellParticleIDsFixed,
                                                         const uint64_t* particleCellIDsSeq,
                                                         const uint*     numParticlesPerCell,
                                                         const uint*     numNeighborsPrefixSums,
                                                         const uint      maxParticlesPerCell,
                                                         const uint      pairListOffset,
                                                         const uint      numObstacles,
                                                         const uint      numParticles,
                                                         const uint      numCells,
                                                         uint2*          pairList)
{
    constexpr uint NUM_NEIGHBOR_CELLS = 27;

    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    // Use sequential buffer: particleCellIDsSeq[tID] = (cellID<<32 | particleID)
    // cellParticleIDsFixed uses 2D layout [cellID*maxPerCell+slot] for neighbour traversal.
    const uint64_t packed = particleCellIDsSeq[tID];
    const uint     cell   = (uint)(packed >> 32);
    const uint     i      = (uint)(packed & 0xFFFFFFFF);

    if(cell == UINT_MAX || cell >= numCells)
        return;

    const uint* neighborCells = &cellNeighborsList[NUM_NEIGHBOR_CELLS * cell];
    uint        insertIndex   = pairListOffset + numNeighborsPrefixSums[i - numObstacles];

    // Loop over all neighboring cells
    for(uint cID = 0; cID < NUM_NEIGHBOR_CELLS; ++cID)
    {
        uint c = neighborCells[cID];

        if(c == UINT_MAX || c < cell)
            continue;

        // Compute cellStart for fixed layout: cellID * maxParticlesPerCell
        uint cellStart = c * maxParticlesPerCell;
        uint cellEnd   = cellStart + numParticlesPerCell[c];

        // Loop through particles in the neighbor cell.
        // For same-cell pairs (c == cell): use i < j to avoid duplicates.
        // For cross-cell pairs (c > cell): the c >= cell condition already ensures
        // only this thread processes the (cell, c) cell pair, so emit all pairs.
        for(uint p = cellStart; p < cellEnd; ++p)
        {
            const uint j = (uint)(cellParticleIDsFixed[p] & 0xFFFFFFFF);

            if(c == cell && i >= j)
                continue;

            // Always write (min,max) to keep pairs canonically ordered.
            pairList[insertIndex++] = (i < j) ? make_uint2(i, j) : make_uint2(j, i);
        }
    }
}