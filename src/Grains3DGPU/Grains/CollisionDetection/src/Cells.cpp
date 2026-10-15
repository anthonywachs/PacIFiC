#include "Cells.hh"
#include "Transform3.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ Cells<T, OrderingScheme>::Cells()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with min and max points along with the extent of each cell.
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__
    Cells<T, OrderingScheme>::Cells(const Vector3<T>& min, const Vector3<T>& max, T cellSize)
    : m_minCorner(min)
    , m_maxCorner(max)
{
    resize(cellSize);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ Cells<T, OrderingScheme>::~Cells()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the min corner point of the linked cell
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ const Vector3<T>& Cells<T, OrderingScheme>::getMinCorner() const
{
    return (m_minCorner);
}

// -------------------------------------------------------------------------------------------------
// Gets the max corner point of the linked cell
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ const Vector3<T>& Cells<T, OrderingScheme>::getMaxCorner() const
{
    return (m_maxCorner);
}

// -------------------------------------------------------------------------------------------------
// Gets the min corner point of the linked cell
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ const Vector3<T>& Cells<T, OrderingScheme>::getMinCornerLinkedCell() const
{
    return (m_minCornerLinkedCell);
}

// -------------------------------------------------------------------------------------------------
// Gets the extent of each cell
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ T Cells<T, OrderingScheme>::getCellSize() const
{
    return (m_cellSize);
}

// -------------------------------------------------------------------------------------------------
// Gets the number of cells along each direction and total
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint4 Cells<T, OrderingScheme>::getNumCellsPerDirection() const
{
    return (m_numCells);
}

// -------------------------------------------------------------------------------------------------
// Gets the number of cells
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::getNumCells() const
{
    return (m_numCells.w);
}

// -------------------------------------------------------------------------------------------------
// Gets the required size for neighbor cells
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::getSizeOfNeighborCells() const
{
    return (27 * m_numCells.w);
}

// -------------------------------------------------------------------------------------------------
// Resizes the linked cells
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ void Cells<T, OrderingScheme>::resize(const T cellSize)
{
    T DX = m_maxCorner[X] - m_minCorner[X];
    T DY = m_maxCorner[Y] - m_minCorner[Y];
    T DZ = m_maxCorner[Z] - m_minCorner[Z];

    m_cellSize     = cellSize;
    m_cellSize_inv = T(1) / cellSize;
    // Use ceiling to ensure the grid is large enough to contain the entire
    // domain
    m_numCells.x = max(1u, uint(ceil(DX * m_cellSize_inv)));
    m_numCells.y = max(1u, uint(ceil(DY * m_cellSize_inv)));
    m_numCells.z = max(1u, uint(ceil(DZ * m_cellSize_inv)));

    // Apply Morton code constraints if using Morton ordering
    if constexpr(OrderingScheme == CellOrdering::MORTON)
    {
        // Morton code limitation: each dimension must be <= 1024 (10 bits)
        constexpr uint MAX_MORTON_DIM = 1024;
        if(m_numCells.x > MAX_MORTON_DIM || m_numCells.y > MAX_MORTON_DIM
           || m_numCells.z > MAX_MORTON_DIM)
        {
            // Scale down to fit Morton code constraints
            T scale = max({T(m_numCells.x), T(m_numCells.y), T(m_numCells.z)}) / T(MAX_MORTON_DIM);
            m_cellSize *= scale;
            m_cellSize_inv = T(1) / m_cellSize;
            m_numCells.x   = max(1u, uint(ceil(DX * m_cellSize_inv)));
            m_numCells.y   = max(1u, uint(ceil(DY * m_cellSize_inv)));
            m_numCells.z   = max(1u, uint(ceil(DZ * m_cellSize_inv)));
        }
    }

    m_numCells.w = m_numCells.x * m_numCells.y * m_numCells.z;

    // Center the domain within the cell grid
    m_minCornerLinkedCell[X] = m_minCorner[X] - (m_numCells.x * m_cellSize - DX) * T(0.5);
    m_minCornerLinkedCell[Y] = m_minCorner[Y] - (m_numCells.y * m_cellSize - DY) * T(0.5);
    m_minCornerLinkedCell[Z] = m_minCorner[Z] - (m_numCells.z * m_cellSize - DZ) * T(0.5);
}

// -------------------------------------------------------------------------------------------------
// Generates neighbor list for cells with ordering-specific implementation
// TODO: We can also improve this so each thread can take more than one cell,
// but that would only be useful for cases where this is a bottleneck.
// TODO: The length of the array is fixed as 27 * m_numCells.w. We can implement
// a more dynamic approach where we first compute the maximum possible
// number of neighbors for each cell, and then use that to allocate the neighbor
// list. This would be more efficient in terms of memory usage.
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ void
    Cells<T, OrderingScheme>::generateNeighborCells(uint* neighborCells, uint start, uint end) const
{
    if(end == 0)
        end = m_numCells.w;
    constexpr uint numNeighbors = 27;

    // Precompute neighbors for each cell
    // clang-format off
    for(uint cellIndex = start; cellIndex < end; ++cellIndex)
    {
        uint offset = numNeighbors * cellIndex;
        // Always use linear decode for row indexing to ensure dense [0..N)
        uint3 cellId = decodeLinearHash(cellIndex);
        
        for(int k = -1; k < 2; ++k) {
        for(int j = -1; j < 2; ++j) {
        for(int i = -1; i < 2; ++i) {
            int nx = cellId.x + i;
            int ny = cellId.y + j;
            int nz = cellId.z + k;

            // Check if the neighbor is within bounds
            if(nx >= 0 && nx < m_numCells.x && 
               ny >= 0 && ny < m_numCells.y && 
               nz >= 0 && nz < m_numCells.z)
            {
                // Store dense linear index for neighbor
                uint neighborIdx = computeLinearHash((uint)nx, (uint)ny, (uint)nz);
                neighborCells[offset++] = neighborIdx;
            }
            else
            {
                // Invalid neighbor
                neighborCells[offset++] = UINT_MAX;
            }
        } } }
    }
    // clang-format on
}

// -------------------------------------------------------------------------------------------------
// Checks if a cell Id is in range
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ bool Cells<T, OrderingScheme>::isValid(const uint3& id) const
{
    return (id.x < m_numCells.x && id.y < m_numCells.y && id.z < m_numCells.z);
}

// -------------------------------------------------------------------------------------------------
// Returns the 3d Id of the cell which the point belongs to
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint3 Cells<T, OrderingScheme>::computeCellID(const Vector3<T>& p,
                                                             bool              checkIfValid) const
{
    const T* __RESTRICT__ pt    = p.getBuffer();
    const T* __RESTRICT__ minPt = m_minCornerLinkedCell.getBuffer();
    uint3                 cellId;
    // static_cast is faster than floor, though it comes with a cost ...
    // if the operand is negative, cast truncates toward zero, while floor
    // rounds down. However, the operand should be always non-negative in our
    // case.
    cellId.x = static_cast<int>((pt[X] - minPt[X]) * m_cellSize_inv);
    cellId.y = static_cast<int>((pt[Y] - minPt[Y]) * m_cellSize_inv);
    cellId.z = static_cast<int>((pt[Z] - minPt[Z]) * m_cellSize_inv);
    // cellId.x = floor((p[X] - m_minCornerLinkedCell[X]) * m_cellSize_inv);
    // cellId.y = floor((p[Y] - m_minCornerLinkedCell[Y]) * m_cellSize_inv);
    // cellId.z = floor((p[Z] - m_minCornerLinkedCell[Z]) * m_cellSize_inv);
    if(checkIfValid)
        GAssert(isValid(cellId), "Linked cell range exceeded!");

    return (cellId);
}

// -------------------------------------------------------------------------------------------------
// Returns the 3d Id of the cell given its hash (ordering-specific)
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint3 Cells<T, OrderingScheme>::computeCellID(const uint cellHash) const
{
    if constexpr(OrderingScheme == CellOrdering::MORTON)
    {
        return decodeMortonCode(cellHash);
    }
    else  // LINEAR ordering
    {
        return decodeLinearHash(cellHash);
    }
}

// -------------------------------------------------------------------------------------------------
// Returns the cell hash value of a given point
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeCellHash(const Vector3<T>& p) const
{
    return (computeCellHash(computeCellID(p)));
}

// -------------------------------------------------------------------------------------------------
// Returns the cell hash value from the 3d Id of the cell (ordering-specific)
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeCellHash(const uint3& cellID) const
{
    if constexpr(OrderingScheme == CellOrdering::MORTON)
    {
        return computeMortonCode(cellID.x, cellID.y, cellID.z);
    }
    else  // LINEAR ordering
    {
        return computeLinearHash(cellID.x, cellID.y, cellID.z);
    }
}

// -------------------------------------------------------------------------------------------------
// Returns the cell hash value from the Id along each axis (ordering-specific)
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeCellHash(const uint i,
                                                              const uint j,
                                                              const uint k) const
{
    if constexpr(OrderingScheme == CellOrdering::MORTON)
    {
        return computeMortonCode(i, j, k);
    }
    else  // LINEAR ordering
    {
        return computeLinearHash(i, j, k);
    }
}

// -------------------------------------------------------------------------------------------------
// Returns dense linear index from point
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeDenseIndex(const Vector3<T>& p,
                                                                bool checkIfValid) const
{
    const uint3 id = computeCellID(p, checkIfValid);
    if(!isValid(id))
        return UINT_MAX;
    return computeLinearHash(id.x, id.y, id.z);
}

// -------------------------------------------------------------------------------------------------
// Returns dense linear index from 3D cell id
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeDenseIndex(const uint3& cellId) const
{
    if(!isValid(cellId))
        return UINT_MAX;
    return computeLinearHash(cellId.x, cellId.y, cellId.z);
}

// -------------------------------------------------------------------------------------------------
// Returns dense linear index from i,j,k
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeDenseIndex(uint i, uint j, uint k) const
{
    if(!(i < m_numCells.x && j < m_numCells.y && k < m_numCells.z))
        return UINT_MAX;
    return computeLinearHash(i, j, k);
}

// -------------------------------------------------------------------------------------------------
// Returns Morton key from dense linear index
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::mortonKeyFromLinearIndex(uint linearIndex) const
{
    // Decode linear to (x,y,z), then compute Morton code
    const uint3 id = decodeLinearHash(linearIndex);
    return computeMortonCode(id.x, id.y, id.z);
}

// -------------------------------------------------------------------------------------------------
// Computes linear hash for 3D coordinates
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeLinearHash(uint x, uint y, uint z) const
{
    return ((z * m_numCells.y + y) * m_numCells.x + x);
}

// -------------------------------------------------------------------------------------------------
// Decodes linear hash to 3D coordinates
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint3 Cells<T, OrderingScheme>::decodeLinearHash(uint hash) const
{
    uint z = hash / (m_numCells.x * m_numCells.y);
    uint y = (hash / m_numCells.x) % m_numCells.y;
    uint x = hash % m_numCells.x;
    return uint3{x, y, z};
}

// -------------------------------------------------------------------------------------------------
// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::expandBits(uint v) const
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// -------------------------------------------------------------------------------------------------
// Compresses a 30-bit Morton-encoded value back to 10 bits
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::compactBits(uint v) const
{
    v &= 0x49249249u;
    v = (v | (v >> 2)) & 0xC30C30C3u;
    v = (v | (v >> 4)) & 0x0F00F00Fu;
    v = (v | (v >> 8)) & 0xFF0000FFu;
    v = (v | (v >> 16)) & 0x000003FFu;
    return v;
}

// -------------------------------------------------------------------------------------------------
// Computes Morton code for 3D coordinates
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint Cells<T, OrderingScheme>::computeMortonCode(uint x, uint y, uint z) const
{
    // Ensure coordinates are within 10-bit range (0-1023)
    x = min(x, 1023u);
    y = min(y, 1023u);
    z = min(z, 1023u);

    return expandBits(x) | (expandBits(y) << 1) | (expandBits(z) << 2);
}

// -------------------------------------------------------------------------------------------------
// Decodes Morton code to 3D coordinates
template <typename T, CellOrdering OrderingScheme>
__HOSTDEVICE__ uint3 Cells<T, OrderingScheme>::decodeMortonCode(uint code) const
{
    uint3 result;
    result.x = compactBits(code);
    result.y = compactBits(code >> 1);
    result.z = compactBits(code >> 2);
    return result;
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation for both ordering schemes
template class Cells<float, CellOrdering::LINEAR>;
template class Cells<double, CellOrdering::LINEAR>;
template class Cells<float, CellOrdering::MORTON>;
template class Cells<double, CellOrdering::MORTON>;
