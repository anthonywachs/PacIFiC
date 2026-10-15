#ifndef _NEIGHBORLIST_HH_
#define _NEIGHBORLIST_HH_

#include "GrainsMemBuffer.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The class NeighborList.

    This class provides functionalities to create a neighbor list for components in the simulation.
    It is used to limit the collision detection to only neighboring components. It is one of the
    main differences between the Grains3D and its GPU version as it is useful to avoid thread
    divergence. This is the base class and derived classes should implement the methods. This design
    gives the flexibility to use different types of neighbor lists. For instance, the neighbor list
    can be created using an O(n^2) algorithm for systems with a small number of components or using
    a more sophisticated algorithm for larger systems.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class NeighborList
{
    static_assert(M == MemType::HOST || M == MemType::DEVICE,
                  "NeighborList only supports MemType::HOST or MemType::DEVICE");

protected:
    /** @name Parameters */
    //@{
    /** \brief Pair list */
    GrainsMemBuffer<uint2, M> m_pairList;
    /** \brief Pair count */
    uint* m_pairCount;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    NeighborList()
    {
        if constexpr(M == MemType::DEVICE)
        {
            cudaErrCheck(cudaMallocManaged(&m_pairCount, sizeof(uint)));
        }
        else
        {
            m_pairCount = new uint;
        }
        *m_pairCount = 0;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    virtual ~NeighborList()
    {
        if constexpr(M == MemType::DEVICE)
        {
            cudaErrCheck(cudaFree(m_pairCount));
        }
        else
        {
            delete m_pairCount;
        }
    }
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Gets pair list */
    const GrainsMemBuffer<uint2, M>& getBuffer() const
    {
        return m_pairList;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets pair list data */
    uint2* getData()
    {
        return m_pairList.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets size of pair list */
    uint getSize() const
    {
        return *m_pairCount;
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the neighbor list
    @param positions array of positions */
    virtual bool updateNeighborList(GrainsMemBuffer<Vector3<T>, M>& positions,
                                    const uint                      nObstacles,
                                    const uint                      nParticles)
        = 0;
    //@}
};

#endif