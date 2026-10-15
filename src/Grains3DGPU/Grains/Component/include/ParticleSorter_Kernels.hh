#ifndef _PARTICLESORTER_KERNELS_HH_
#define _PARTICLESORTER_KERNELS_HH_

#include "Basic.hh"
#include "Cells.hh"

// =================================================================================================
/** @brief The class ParticleSorter_Kernels.

    This header file contains the declarations of the various kernels used for
    sorting particles based on Morton codes.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name ParticleSorter_Kernels: External Kernels */
//@{
/** @brief Deletes the device-heap-allocated Cells object at index 0.
    The Cells object must have been created with device-side 'new' (e.g., via createCells_Kernel).
    Host-side cudaFree cannot free device-heap allocations, so this kernel is required.
    @param cells device pointer array; cells[0] is deleted and set to nullptr */
template <typename T, CellOrdering CO>
__GLOBAL__ void deleteCells_Kernel(Cells<T, CO>** cells)
{
    if(blockIdx.x == 0 && threadIdx.x == 0 && cells[0] != nullptr)
    {
        delete cells[0];
        cells[0] = nullptr;
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes Morton codes for all particles
    @param cells pointer to the Cells object with Morton ordering
    @param positions particle positions
    @param mortonCodes output Morton codes
    @param numParticles number of particles */
template <typename T>
__GLOBAL__ void computeMortonCodes_Kernel(Cells<T, CellOrdering::MORTON>** cells,
                                          const Vector3<T>*                positions,
                                          uint64_t*                        mortonCodes,
                                          uint                             numParticles)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    mortonCodes[tID] = (*cells)->computeCellHash(positions[tID]);
}

// -------------------------------------------------------------------------------------------------
/** @brief Gathers data according to sorted indices (templated)
    @param src source array
    @param dst destination array
    @param indices sorted indices
    @param numParticles number of particles */
template <typename T>
__GLOBAL__ void gather_Kernel(const T* src, T* dst, const uint* indices, uint numParticles)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID >= numParticles)
        return;

    dst[tID] = src[indices[tID]];
}
//@}

#endif
