#ifndef _NEIGHBORLIST_KERNELS_HH_
#define _NEIGHBORLIST_KERNELS_HH_

#include "GrainsMemBuffer.hh"

// =================================================================================================
/** @brief The class NeighborList_Kernels.

    This header file contains the declarations of the various kernels used for
    updating the neighbor list in the simulation. These kernels support both
    host and device.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name NeighborList_Kernels: External Kernels */
//@{
/** @brief Updates the neighbor list on host using an O(n^2) algorithm.
    @param nObstacles number of obstacles
    @param nParticles number of particles
    @param pairList array of pairs */
__HOST__ void
    updateNeighborList_Nsq_Host(const uint nObstacles, const uint nParticles, uint2* pairList);

/** @brief Updates the neighbor list on device using an O(n^2) algorithm.
    @param nObstacles number of obstacles
    @param nParticles number of particles
    @param pairList array of pairs */
__GLOBAL__ void
    updateNeighborList_Nsq_Device(const uint nObstacles, const uint nParticles, uint2* pairList);

/** @brief Updates the neighbor list on host using a linked cell approach.
    @param cellNeighborsList array of neighboring cells for each cell
    @param obstacleIDs array of obstacle IDs
    @param obstacleCellIDs array of obstacle cell IDs
    @param particleIDs array of particle IDs
    @param cellIDs array of cell IDs
    @param cellParticles vector of lists containing particle IDs for each cell
    @param maxCellsPerObstacle maximum number of cells per obstacle
    @param numObstacles number of obstacles
    @param numParticles number of particles
    @param pairList reference to GrainsMemBuffer that will be resized if needed */
__HOST__ void updateNeighborList_LC_Host(const uint*                            cellNeighborsList,
                                         const uint2*                           obstacleIDs,
                                         const uint*                            obstacleCellIDs,
                                         const uint*                            particleIDs,
                                         const uint*                            cellIDs,
                                         const std::vector<std::list<uint>>&    cellParticles,
                                         const uint                             maxCellsPerObstacle,
                                         const uint                             numObstacles,
                                         const uint                             numParticles,
                                         GrainsMemBuffer<uint2, MemType::HOST>& pairList);

/** @brief Generate obstacle-particle pairs on device.
    @param obstacleIDs array of obstacle IDs and cell counts
    @param obstacleCellIDs array of obstacle cell IDs
    @param cellParticleIDs packed uint64 array (u32 bits = cellID, l32 bits = particleID)
    @param cellPrefixSums prefix sum array for cell particle positions
    @param maxCellsPerObstacle maximum number of cells per obstacle
    @param numObstacles number of obstacles
    @param numParticles number of particles
    @param numCells number of cells
    @param pairList array of pairs
    @param pairCount pointer to device memory for storing the total pair count */
__GLOBAL__ void generateObstacleParticlePairs_Device(const uint2*    obstacleIDs,
                                                     const uint*     obstacleCellIDs,
                                                     const uint64_t* cellParticleIDs,
                                                     const uint*     cellPrefixSums,
                                                     const uint      maxCellsPerObstacle,
                                                     const uint      numObstacles,
                                                     const uint      numParticles,
                                                     const uint      numCells,
                                                     uint2*          pairList,
                                                     uint*           pairCount);

/** @brief Count neighbors per particle using linked cells with packed keys.
    @param cellNeighborsList array of neighboring cells for each cell
    @param cellParticleIDs packed uint64 array (upper 32 bits = cellID, lower 32 bits = particleID)
    @param cellPrefixSums array of start IDs for each cell
    @param numParticles number of particles
    @param numObstacles number of obstacles (subtracted from global ID to get local index)
    @param numCells number of cells
    @param neighborCounts output array of neighbor counts per particle */
__GLOBAL__ void countNeighbors_Device(const uint*     cellNeighborsList,
                                      const uint64_t* cellParticleIDs,
                                      const uint*     cellPrefixSums,
                                      const uint      numParticles,
                                      const uint      numObstacles,
                                      const uint      numCells,
                                      uint*           neighborCounts);

/** @brief Updates the neighbor list on device keys.
    @param cellNeighborsList array of neighboring cells for each cell
    @param cellParticleIDs packed uint64 array (upper 32 bits = cellID, lower 32 bits = particleID)
    @param cellPrefixSums array of start IDs for each cell
    @param numNeighborsPrefixSums array of prefix sums of neighbors
    @param pairListOffset offset in pairList where particle-particle pairs should start
    @param numObstacles number of obstacles
    @param numParticles number of particles
    @param numCells number of cells
    @param pairList array of pairs */
__GLOBAL__ void updateNeighborList_LC_Device(const uint*     cellNeighborsList,
                                             const uint64_t* cellParticleIDs,
                                             const uint*     cellPrefixSums,
                                             const uint*     numNeighborsPrefixSums,
                                             const uint      pairListOffset,
                                             const uint      numObstacles,
                                             const uint      numParticles,
                                             const uint      numCells,
                                             uint2*          pairList);

/** @brief Count neighbors per particle for AtomicFixed (direct count lookup).
    @param cellNeighborsList array of neighboring cells for each cell
    @param cellParticleIDs packed uint64 array (upper 32 bits = cellID, lower 32 bits = particleID)
    @param numParticlesPerCell array of particle counts per cell
    @param maxParticlesPerCell maximum particles per cell
    @param numParticles number of particles
    @param numObstacles number of obstacles (subtracted from global ID to get local index)
    @param numCells number of cells
    @param neighborCounts output array of neighbor counts per particle */
__GLOBAL__ void countNeighbors_AtomicFixed_Device(const uint*     cellNeighborsList,
                                                  const uint64_t* cellParticleIDsFixed,
                                                  const uint64_t* particleCellIDsSeq,
                                                  const uint*     numParticlesPerCell,
                                                  const uint      maxParticlesPerCell,
                                                  const uint      numParticles,
                                                  const uint      numObstacles,
                                                  const uint      numCells,
                                                  uint*           neighborCounts);

/** @brief Updates the neighbor list for AtomicFixed (direct count lookup).
    @param cellNeighborsList array of neighboring cells for each cell
    @param cellParticleIDs packed uint64 array (upper 32 bits = cellID, lower 32 bits = particleID)
    @param numParticlesPerCell array of particle counts per cell
    @param numNeighborsPrefixSums array of prefix sums of neighbors
    @param maxParticlesPerCell maximum particles per cell
    @param pairListOffset offset in pairList where particle-particle pairs should start
    @param numObstacles number of obstacles
    @param numParticles number of particles
    @param numCells number of cells
    @param pairList array of pairs */
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
                                                         uint2*          pairList);
//@}

#endif