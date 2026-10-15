#ifndef _CELLSFACTORY_HH_
#define _CELLSFACTORY_HH_

#include <cuda_runtime.h>

#include "Cells.hh"
#include "GrainsMemBuffer.hh"

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
/** @brief GPU kernel to construct the Cells on device. This is mandatory as we cannot access
    device memory addresses on the host. So, we pass a device memory address to a kernel. Memory
    address is then populated within the kernel. */
template <typename T, CellOrdering CO, typename... Arguments>
__GLOBAL__ void createCells_Kernel(Cells<T, CO>** cells,
                                   uint           index,
                                   T              minX,
                                   T              minY,
                                   T              minZ,
                                   T              maxX,
                                   T              maxY,
                                   T              maxZ,
                                   T              size,
                                   uint*          numCells)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID > 0)
        return;

    cells[index]
        = new Cells<T, CO>(Vector3<T>(minX, minY, minZ), Vector3<T>(maxX, maxY, maxZ), size);
    *numCells = cells[index]->getNumCells();
}

// =================================================================================================
/** @brief The class CellsFactory.

    Creates cells object for the simulation.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T, CellOrdering CO = CellOrdering::LINEAR>
class CellsFactory
{
private:
    /**@name Contructors & Destructor */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor (forbidden) */
    CellsFactory() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor (forbidden) */
    ~CellsFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Creates and returns a buffer of cells with memory type handling
        @param minCorner Minimum corner of the domain
        @param maxCorner Maximum corner of the domain
        @param cellSize Size of each cell
        @param cells Memory buffer for storing the cell object(s) */
    template <MemType M>
    static uint create(const Vector3<T>&                  minCorner,
                       const Vector3<T>&                  maxCorner,
                       const T                            cellSize,
                       GrainsMemBuffer<Cells<T, CO>*, M>& cells)
    {
        cells.reserve(1);
        if constexpr(M == MemType::HOST)
        {
            return create_host(minCorner, maxCorner, cellSize, cells);
        }
        else if constexpr(M == MemType::DEVICE)
        {
            return create_device(minCorner, maxCorner, cellSize, cells);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Creates and returns a buffer of cells on host
        @param minCorner Minimum corner of the domain
        @param maxCorner Maximum corner of the domain
        @param cellSize Size of each cell
        @param cells Memory buffer for storing the cell object(s) */
    static uint create_host(const Vector3<T>&                              minCorner,
                            const Vector3<T>&                              maxCorner,
                            const T                                        cellSize,
                            GrainsMemBuffer<Cells<T, CO>*, MemType::HOST>& cells)
    {
        // Safety check
        GAssert(cellSize > 0, "Cell size must be positive! Aborting Grains!");
        cells.initialize(1);
        cells[0] = new Cells<T, CO>(minCorner, maxCorner, cellSize);
        return cells[0]->getNumCells();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Creates and returns a buffer of cells on device
        @param minCorner Minimum corner of the domain
        @param maxCorner Maximum corner of the domain
        @param cellSize Size of each cell
        @param cells Memory buffer for storing the cell object(s) */
    static uint create_device(const Vector3<T>&                                minCorner,
                              const Vector3<T>&                                maxCorner,
                              const T                                          cellSize,
                              GrainsMemBuffer<Cells<T, CO>*, MemType::DEVICE>& cells)
    {
        GrainsMemBuffer<Cells<T, CO>*, MemType::HOST> h_cells(1);
        uint numCells = create_host(minCorner, maxCorner, cellSize, h_cells);
        copyHostToDevice(h_cells, cells);
        // Free the host buffer
        delete h_cells[0];
        return numCells;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Cells objects must be instantiated on device, if we want to use them on device.
        Copying from host is not supported due to runtime polymorphism for this class. This
        function reads a host-side Cells object, and mimics it in a given device buffer. It calls a
        device kernel that is implemented in the source file.
        @param h_LC Host-side Cells object
        @param d_LC Device-side Cells object */
    static void copyHostToDevice(GrainsMemBuffer<Cells<T, CO>*, MemType::HOST>&   h_LC,
                                 GrainsMemBuffer<Cells<T, CO>*, MemType::DEVICE>& d_LC)
    {
        // Allocate the device memory for the linked cells
        d_LC.initialize(h_LC.getSize());
        uint  h_numCells = 0;
        uint* d_numCells;
        cudaMalloc(&d_numCells, sizeof(uint));
        for(uint i = 0; i < h_LC.getSize(); ++i)
        {
            if(h_LC[i] == nullptr)
                continue;

            // Extracting info from the host side object
            Vector3<T> origin        = h_LC[i]->getMinCorner();
            Vector3<T> maxCoordinate = h_LC[i]->getMaxCorner();
            T          size          = h_LC[i]->getCellSize();
            // Safety check
            GAssert(size > 0, "Cell size must be positive! Aborting Grains!");
            createCells_Kernel<<<1, 1>>>(d_LC.getData(),
                                         i,
                                         origin[X],
                                         origin[Y],
                                         origin[Z],
                                         maxCoordinate[X],
                                         maxCoordinate[Y],
                                         maxCoordinate[Z],
                                         size,
                                         d_numCells);
            cudaErrCheck(cudaMemcpy(&h_numCells, d_numCells, sizeof(uint), cudaMemcpyDeviceToHost));
            // GoutWI(9, "LinkedCell with", h_numCells, "cells is created on device.");
        }
        cudaErrCheck(cudaDeviceSynchronize());
        cudaErrCheck(cudaFree(d_numCells));
    }
    //@}
};

#endif
