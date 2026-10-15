#ifndef _GRAINSPARAMETERS_HH_
#define _GRAINSPARAMETERS_HH_

#include "GrainsTimer.hh"
#include "Vector3.hh"

/** @name Enumerations */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Type of neighbor list */
enum class NeighborListType
{
    /** @brief N-Squared neighbor list */
    NSQ = 0,
    /** @brief Linked cell neighbor list */
    LINKEDCELL = 1
};

// -------------------------------------------------------------------------------------------------
/** @brief Type of linked cell */
enum class LinkedCellType
{
    /** @brief Host linked cells */
    HOST = 0,
    /** @brief Sort-based linked cells for device */
    SORTBASED = 1,
    /** @brief Atomic linked cells for device */
    ATOMIC = 2,
    /** @brief Atomic fixed-size linked cells for device (pre-sized 2D array) */
    ATOMICFIXED = 3
};

// -------------------------------------------------------------------------------------------------
/** @brief Type of bounding volume */
enum class BoundingVolumeType
{
    /** @brief No bounding volume */
    OFF = 0,
    /** @brief Oriented Bounding Box */
    OBB = 1,
    /** @brief Oriented Bounding Cylinder */
    OBC = 2
};

// -------------------------------------------------------------------------------------------------
/** @brief Type of narrow-phase detection */
enum class NarrowPhaseType
{
    /** @brief GJK with Johnson's sub-algorithm */
    GJK = 0,
    /** @brief GJK with Signed Volume sub-algorithm */
    GJK_SV = 1
};
//@}

/** @name Structs */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Plain-old-data struct grouping all simulation-level component counts. */
struct ComponentCounts
{
    uint numObstacles  = 0;  ///< Number of obstacles
    uint numParticles  = 0;  ///< Number of moving particles
    uint numPairs      = 0;  ///< Number of active contact pairs (updated each step)
    uint numComposites = 0;  ///< Number of composite bodies
    uint numSubBodies  = 0;  ///< Total sub-body slots across all composites
};

// -------------------------------------------------------------------------------------------------
/** @brief Parameters to track dynamic simulation state during runtime. */
template <typename T>
struct SimulationState
{
    /** @brief Current simulation time. */
    T time = 0;
    /** @brief Number of times neighbor list has been updated. */
    uint neighborListUpdateCount = 0;
    /** @brief did obstacles move in the last step? */
    bool obstaclesMoved = true;
    /** @brief did particles get sorted in the last step? */
    bool particlesSorted = true;
};

// -------------------------------------------------------------------------------------------------
/** @brief Parameters for linked cell configuration */
template <typename T>
struct LinkedCellParameters
{
    /** \brief Minimum corner of the linked cell domain. */
    Vector3<T> minCorner = Vector3<T>(0, 0, 0);
    /** \brief Maximum corner of the linked cell domain. */
    Vector3<T> maxCorner = Vector3<T>(0, 0, 0);
    /** \brief Type of linked cell. */
    LinkedCellType type = LinkedCellType::HOST;
    /** \brief Minimum linked cell size. */
    T minCellSize = 0;
    /** \brief Linked cell size factor. */
    T cellSizeFactor = 1;
    /** \brief Maximum circumscribed radius among all obstacles (used for LC sizing). */
    T maxObstacleRadius = 0;
    /** \brief Maximum number of cells that one obstacle can occupy. */
    uint maxNumCellsPerObstacle = 0;
    /** \brief Maximum number of particles per cell (used for ATOMICFIXED LinkedCell type). */
    uint maxParticlesPerCell = 64;
    /** \brief Initial number of pairs per particle. */
    uint initialNumberOfPairsPerParticle = 16;
    /** \brief If using adaptive skin, this is the desired number of iterations that the skin should
        be valid for. If it is set to 0, then we don't use adaptive skin. If set to 1, the skin is
        updated every iteration and trivially it gives worse performance. */
    uint updateFrequency = 0;
    /** \brief Number of iterations between each sorting. */
    uint sortFrequency = 0;
};

// -------------------------------------------------------------------------------------------------
/** @brief Parameters for collision detection configuration. */
template <typename T>
struct CollisionDetectionParameters
{
    /** \brief Type of neighbor list. */
    NeighborListType neighborListType = NeighborListType::NSQ;
    /** \brief LinkedCell parameters. */
    LinkedCellParameters<T> linkedCellParameters;
    /** \brief Type of bounding volume. */
    BoundingVolumeType boundingVolumeType = BoundingVolumeType::OFF;
    /** \brief Type of narrow-phase detection. */
    NarrowPhaseType narrowPhaseType = NarrowPhaseType::GJK;
    /** \brief Use GJK warm-start acceleration. */
    bool gjkAcceleration = false;
    /** \brief If true, relative transformations (B in A-local frame) are computed once before
        narrow-phase detection and cached; set to false (e.g. for spheres) to skip that step and
        call the GJK overload that works directly in world frame instead. */
    bool useRelativeTransformations = true;
    /** \brief Use pre-built ShapeData array for vtable-free GJK support evaluation. */
    bool usePrebuiltShapes = false;
};
//@}

// =================================================================================================
/** @brief Global Parameters in Grains.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class GrainsParameters
{
public:
    /** @name Parameters */
    //@{
    /* GPU */
    /** \brief is simulation on GPU? */
    static bool m_isGPU;
    /** \brief GPU device properties */
    static cudaDeviceProp m_GPU;

    /* Dynamic Settings */
    /** \brief Simulation state */
    static SimulationState<T> m_simulationState;

    /* Spatial */
    /** @brief Global domain origin */
    static Vector3<T> m_origin;
    /** @brief Global domain dimension */
    static Vector3<T> m_maxCoordinate;
    /** @brief Is simulation periodic? */
    static bool m_isPeriodic;

    /* Collision Detection */
    /** \brief Collision detection parameters. */
    static CollisionDetectionParameters<T> m_collisionDetection;

    /* Material */
    /** \brief Map from material name to an uint ID */
    static std::unordered_map<std::string, uint> m_materialMap;
    /** \brief Number of different possible contact pairs (incl. obs-obs) */
    static uint m_numContactPairs;
    /** \brief Is contact with memory activated? */
    static bool m_isContactWithMemory;
    /** \brief Use compaction to work only on active pairs on GPU.
        When false the force kernel runs over all pairs. */
    static bool m_useCompaction;

    /* Temporal */
    /** @brief Initial simulation time */
    static T m_tStart;
    /** @brief End simulation time */
    static T m_tEnd;
    /** @brief Simulation time step */
    static T m_dt;
    /** @brief is time integrator leap-frog? */
    static bool m_isLeapFrog;

    /* Physical */
    /** \brief Gravity vector */
    static Vector3<T> m_gravity;

    /* Post-Processing */
    /** \brief Queue of simulation time to write Post-Processing */
    static std::queue<T> m_tSave;
    /** @brief Frequency of time output (print every N steps, 0=never, 1=every step) */
    static uint m_verbosityFrequency;
    /** \brief Simulation-loop timer */
    static GrainsSimTimer m_simTimer;
    /** \brief CDM sub-stage timer */
    static GrainsCDMTimer m_cdmTimer;
    /** \brief ForceModule sub-stage timer */
    static GrainsForceTimer m_fmTimer;
    //@}
};

#endif