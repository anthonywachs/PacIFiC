#ifndef _NEIGHBORLISTFACTORY_HH_
#define _NEIGHBORLISTFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "NeighborList.hh"
#include "NeighborList_LinkedCell.hh"
#include "NeighborList_Nsq.hh"

// =================================================================================================
/** @brief The class NeighborListFactory.

    Creates the neighbor list for the simulation.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class NeighborListFactory
{
    static_assert(M == MemType::HOST || M == MemType::DEVICE,
                  "NeighborListFactory only supports MemType::HOST or MemType::DEVICE");

private:
    /**@name Contructors & Destructor */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor (forbidden) */
    NeighborListFactory() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor (forbidden) */
    ~NeighborListFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Creates and returns a unique_ptr to a NeighborList object
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param CD Collision detection parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles
        @return unique_ptr to the created neighbor list */
    static std::unique_ptr<NeighborList<T, M>>
        create(const GrainsMemBuffer<RigidBody<T>*, M>* rb,
               const GrainsMemBuffer<Vector3<T>, M>&    positions,
               const GrainsMemBuffer<Quaternion<T>, M>& quaternions,
               const CollisionDetectionParameters<T>&   CD,
               const uint                               nObstacles,
               const uint                               nParticles)
    {
        // Assertions
        GAssert(rb->getSize() == nObstacles + nParticles, "Rigid body size mismatch");
        GAssert(positions.getSize() == nObstacles + nParticles, "Positions size mismatch");
        GAssert(quaternions.getSize() == nObstacles + nParticles, "Quaternions size mismatch");

        // Global parameters
        NeighborListType type = CD.neighborListType;

        std::unique_ptr<NeighborList<T, M>> NL;
        if(type == NeighborListType::NSQ)
        {
            NL = std::make_unique<NeighborList_Nsq<T, M>>(nObstacles, nParticles);
        }
        else if(type == NeighborListType::LINKEDCELL)
        {
            NL = std::make_unique<NeighborList_LinkedCell<T, M>>(rb,
                                                                 positions,
                                                                 quaternions,
                                                                 CD.linkedCellParameters,
                                                                 nObstacles,
                                                                 nParticles);
        }

        // Sanity check to ensure neighbor list was created
        GAssert(NL != nullptr, "Neighbor list creation failed.");
        return NL;
    }
    //@}
};

#endif
