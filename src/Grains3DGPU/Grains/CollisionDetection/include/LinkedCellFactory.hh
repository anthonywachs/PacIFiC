#ifndef _LINKEDCELLFACTORY_HH_
#define _LINKEDCELLFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "LinkedCell.hh"
#include "LinkedCell_Atomic.hh"
#include "LinkedCell_AtomicFixed.hh"
#include "LinkedCell_Host.hh"
#include "LinkedCell_SortBased.hh"

// =================================================================================================
/** @brief The class LinkedCellFactory.

    Creates the linked cell structure for the simulation.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T, MemType M>
class LinkedCellFactory
{
    static_assert(M == MemType::HOST || M == MemType::DEVICE,
                  "LinkedCellFactory only supports MemType::HOST or MemType::DEVICE");

private:
    /** @name Constructors & Destructor */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor (forbidden) */
    LinkedCellFactory() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor (forbidden) */
    ~LinkedCellFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Creates and returns a buffer of LinkedCell objects
        @param rb Rigid body buffer
        @param positions Positions buffer
        @param quaternions Quaternions buffer
        @param linkedCellParameters Linked cell parameters
        @param nObstacles number of obstacles
        @param nParticles number of particles
        @return unique_ptr to the created LinkedCell object */
    static std::unique_ptr<LinkedCell<T, M>>
        create(const GrainsMemBuffer<RigidBody<T>*, M>* rb,
               const GrainsMemBuffer<Vector3<T>, M>&    positions,
               const GrainsMemBuffer<Quaternion<T>, M>& quaternions,
               const LinkedCellParameters<T>&           linkedCellParameters,
               const uint                               nObstacles,
               const uint                               nParticles)
    {
        auto type = linkedCellParameters.type;
        // Create the linked cell object
        if constexpr(M == MemType::HOST)
        {
            return std::make_unique<LinkedCell_Host<T>>(rb,
                                                        positions,
                                                        quaternions,
                                                        linkedCellParameters,
                                                        nObstacles,
                                                        nParticles);
        }
        else if constexpr(M == MemType::DEVICE)
        {
            if(type == LinkedCellType::SORTBASED)
            {
                return std::make_unique<LinkedCell_SortBased<T>>(rb,
                                                                 positions,
                                                                 quaternions,
                                                                 linkedCellParameters,
                                                                 nObstacles,
                                                                 nParticles);
            }
            else if(type == LinkedCellType::ATOMIC)
            {
                return std::make_unique<LinkedCell_Atomic<T>>(rb,
                                                              positions,
                                                              quaternions,
                                                              linkedCellParameters,
                                                              nObstacles,
                                                              nParticles);
            }
            else if(type == LinkedCellType::ATOMICFIXED)
            {
                return std::make_unique<LinkedCell_AtomicFixed<T>>(rb,
                                                                   positions,
                                                                   quaternions,
                                                                   linkedCellParameters,
                                                                   nObstacles,
                                                                   nParticles);
            }
            else
                GAbort("LinkedCell type not supported on device. Aborting "
                       "Grains!");
        }

        return nullptr;  // Should never reach here
    }
    //@}
};

#endif
