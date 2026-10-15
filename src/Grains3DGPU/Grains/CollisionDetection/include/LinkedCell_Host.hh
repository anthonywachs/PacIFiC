#ifndef _LINKEDCELL_HOST_HH_
#define _LINKEDCELL_HOST_HH_

#include <cstring>
#include <list>
#include <unordered_map>
#include <vector>

#include "Cells.hh"
#include "CellsFactory.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "LinkedCell.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The class LinkedCell_Host.

    This class provides functionalities to manage linked cells for collision detection in the
    simulation on the host. It uses std::vector of std::list for each cell to efficiently manage
    particle assignments. When updating, it checks if particle cell IDs have changed and moves
    particles between cells accordingly.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
template <typename T>
class LinkedCell_Host : public LinkedCell<T, MemType::HOST>
{
    using LC = LinkedCell<T, MemType::HOST>;
    using LC::m_cells;
    using LC::m_neighborCells;
    using LC::m_numCells;
    using LC::m_numObstacles;
    using LC::m_numParticles;
    using LC::m_positions;
    using LC::m_useAdaptiveSkin;

private:
    /** @name Host-specific storage */
    //@{
    /** \brief Buffer of particle IDs (Host-specific: needed for iteration) */
    GrainsMemBuffer<uint, MemType::HOST> m_particleID;
    /** \brief Buffer of cells that particles belong to (Host-specific: needed for updates) */
    GrainsMemBuffer<uint, MemType::HOST> m_cellID;
    /** \brief Temporary buffer to store old particle hashes during updates */
    GrainsMemBuffer<uint, MemType::HOST> m_oldCellID;
    /** \brief Vector of lists, one list per cell containing particle IDs */
    std::vector<std::list<uint>> m_cellParticles;
    /** \brief Map to store iterators to particle positions in cell lists for
        O(1) removal */
    std::unordered_map<uint, std::list<uint>::iterator> m_particleIteratorMap;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles */
    LinkedCell_Host(const GrainsMemBuffer<RigidBody<T>*, MemType::HOST>* rb,
                    const GrainsMemBuffer<Vector3<T>, MemType::HOST>&    positions,
                    const GrainsMemBuffer<Quaternion<T>, MemType::HOST>& quaternions,
                    const LinkedCellParameters<T>&                       linkedCellParameters,
                    const uint                                           nObstacles,
                    const uint                                           nParticles)
        : LinkedCell<T, MemType::HOST>(
              rb, positions, quaternions, linkedCellParameters, nObstacles, nParticles)
        , m_particleID(nParticles)
        , m_cellID(nParticles)
        , m_oldCellID(nParticles)
    {
        // Initialize particle and cell ID buffers
        m_particleID.sequence(m_numObstacles);
        m_cellID.fill(UINT_MAX);
        m_oldCellID.fill(UINT_MAX);

        // Initialize vector of lists for each cell
        m_cellParticles.resize(m_numCells);
        // Reserve space in iterator map for efficiency
        m_particleIteratorMap.reserve(nParticles);

        // Populate initial cell assignments
        populateInitialCells();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~LinkedCell_Host() = default;
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Gets particle IDs */
    const uint* getParticleIDs() const override
    {
        return m_particleID.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets cell IDs */
    const uint* getCellIDs() const override
    {
        return m_cellID.getData();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets all cell particle lists */
    const std::vector<std::list<uint>>& getCellParticles() const
    {
        return m_cellParticles;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Gets particles in a specific cell
        @param cellID the cell ID */
    const std::list<uint>& getParticlesInCell(uint cellID) const
    {
        return m_cellParticles[cellID];
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Populates initial cell assignments for all components */
    void populateInitialCells()
    {
        for(uint i = 0; i < m_numParticles; ++i)
        {
            uint particleID = m_particleID[i];
            uint cellID     = m_cellID[i];
            addParticleToCell(particleID, cellID);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Adds a particle to a specific cell
        @param particleID the particle ID
        @param cellID the cell ID */
    void addParticleToCell(uint particleID, uint cellID)
    {
        if(cellID == UINT_MAX)
            return;

        // Add particle to the cell's list
        m_cellParticles[cellID].push_front(particleID);

        // Store the iterator for O(1) removal later
        m_particleIteratorMap[particleID] = m_cellParticles[cellID].begin();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Removes a particle from its current cell
        @param particleID the particle ID
        @param cellID the cell ID */
    void removeParticleFromCurrentCell(uint particleID, uint cellID)
    {
        if(cellID == UINT_MAX)
            return;

        auto iterIt = m_particleIteratorMap.find(particleID);

        if(iterIt != m_particleIteratorMap.end())
        {
            // Remove from the cell's list using the stored iterator
            m_cellParticles[cellID].erase(iterIt->second);

            // Clean up iterator map
            m_particleIteratorMap.erase(iterIt);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Moves a particle from one cell to another
        @param particleID the particle ID
        @param oldCellID the old cell ID
        @param newCellID the new cell ID */
    void moveParticleToCell(uint particleID, uint oldCellID, uint newCellID)
    {
        removeParticleFromCurrentCell(particleID, oldCellID);
        addParticleToCell(particleID, newCellID);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Updates the linked cells based on particle transformations */
    bool updateLinkedCells() override
    {
        auto& SS = GrainsParameters<T>::m_simulationState;

        // Store old cell IDs before updating
        m_oldCellID.copyFrom(m_cellID);

        // Re-link obstacles if they moved
        bool hasChanges = false;
        if(SS.obstaclesMoved)
        {
            this->linkObstacles();
            hasChanges = true;
        }

        // Update particle cell IDs and move particles that changed cells
        const Vector3<T>* p = m_positions->getData() + m_numObstacles;
        for(uint i = 0; i < m_numParticles; ++i)
        {
            uint particleID = m_particleID[i];
            uint oldCellID  = m_oldCellID[i];
            uint newCellID  = m_cells[0]->computeCellHash(p[i]);
            m_cellID[i]     = newCellID;
            if(oldCellID != newCellID)
            {
                moveParticleToCell(particleID, oldCellID, newCellID);
                hasChanges = true;
            }
        }

        return hasChanges;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Collects particle IDs in the candidate's cell and its neighbors.
        Appends IDs less than maxIndex into out.
        @param candidate candidate world-space position
        @param maxIndex only IDs < maxIndex are considered existing
        @param out output buffer to append IDs */
    void collectPotentialNeighbors(const Vector3<T>&  candidate,
                                   uint               maxIndex,
                                   std::vector<uint>& out) const
    {
        constexpr uint  NUM_NEIGHBOR_CELLS = 27;  // Number of neighboring cells
        const Cells<T>* cells              = this->getLinkedCell()[0];
        const uint      candidateCellID    = cells->computeCellHash(candidate);

        // Build the 27-cell neighborhood of the candidate once (used for both obstacles and
        // particles so that the search radius is symmetric in both cases)
        const uint* allNeighbors  = this->getCellNeighborsList();
        const uint* neighborCells = &allNeighbors[NUM_NEIGHBOR_CELLS * candidateCellID];

        // Collect obstacles: check the candidate's own cell AND all 27 neighboring cells
        // against each obstacle's registered cell list.  The particle-particle search below
        // uses the same 27-cell radius, so the two searches are now symmetric -- without this
        // a candidate sitting in a cell adjacent to an obstacle boundary would miss the
        // obstacle entirely and get inserted overlapping it.
        const uint2* obstacleIDs         = this->getObstacleIDs();
        const uint*  obstacleCellIDs     = this->getObstacleCellIDs();
        const uint   maxCellsPerObstacle = this->getMaxCellsPerObstacle();

        for(uint i = 0; i < m_numObstacles; ++i)
        {
            const uint obstacleIndex      = obstacleIDs[i].x;
            const uint numCellsToTraverse = obstacleIDs[i].y;
            const uint offset             = i * maxCellsPerObstacle;

            bool found = false;

            // Check candidate's own cell
            for(uint c = 0; c < numCellsToTraverse && !found; ++c)
            {
                if(obstacleCellIDs[offset + c] == candidateCellID)
                {
                    out.push_back(obstacleIndex);
                    found = true;
                }
            }

            // Check all 27 neighboring cells of the candidate
            for(uint n = 0; n < NUM_NEIGHBOR_CELLS && !found; ++n)
            {
                const uint nc = neighborCells[n];
                if(nc == UINT_MAX || nc == candidateCellID)
                    continue;
                for(uint c = 0; c < numCellsToTraverse && !found; ++c)
                {
                    if(obstacleCellIDs[offset + c] == nc)
                    {
                        out.push_back(obstacleIndex);
                        found = true;
                    }
                }
            }
        }

        // Same-cell particles
        const auto& currentCellParticles = m_cellParticles[candidateCellID];
        if(!currentCellParticles.empty())
        {
            for(const uint p : currentCellParticles)
                if(p < maxIndex)
                    out.push_back(p);
        }

        // Neighbor cells for particles
        for(uint n = 0; n < NUM_NEIGHBOR_CELLS; ++n)
        {
            const uint c = neighborCells[n];
            if(c == UINT_MAX || c == candidateCellID)
                continue;
            const auto& neighborCellParticles = m_cellParticles[c];
            for(const uint p : neighborCellParticles)
                if(p < maxIndex)
                    out.push_back(p);
        }
    }
    //@}
};

#endif