#ifndef _BENCHMARKCONFIG_HH_
#define _BENCHMARKCONFIG_HH_

#include "GrainsParameters.hh"
#include "Vector3.hh"

// -------------------------------------------------------------------------------------------------
/** @brief Floating-point precision for the benchmark run. */
enum class PrecisionType
{
    SINGLE = 0,
    DOUBLE = 1
};

// -------------------------------------------------------------------------------------------------
/** @brief Particle convex shape type. */
enum class ParticleShapeType
{
    SPHERE       = 0,
    BOX          = 1,
    SUPERQUADRIC = 2
};

// =================================================================================================
/** @brief All parameters that fully define one benchmark scenario.
    These are passed to BenchmarkRunner<T, M>, which sets the GrainsParameters static fields,
    constructs a CollisionDetectionModule, warms it up, then measures the per-stage pipeline
    timings using the built-in GrainsParameters<T>::m_cdmTimer counter.

    === Timing modes (numWarmupCalls / numMeasureCalls) ===
      numWarmupCalls = 0   :  The very first run() call (NL rebuild + possible sort) is
                              included in the timing.  Use this to benchmark cold-path cost
                              (NL construction, first-call sort).
      numWarmupCalls >= 1  :  All warmup calls are executed before the timer is reset, so the
                              NL is already stable before timing starts.  This gives steady-state
                              narrow-phase / BV / transform cost per call.

    === BV + relative-transform interaction ===
      The OBB pre-filter (filterPairsBV) is only executed on the relative-transform path.
      Setting BV=OBB with useRelativeTransformations=false is valid but the BV stage is a
      no-op at runtime; bvFilterTime will be 0.

    === Sorting ===
      sortFrequency N > 0  :  sort triggers when neighborListUpdateCount % N == 0.
                              Since counts only increment on NL rebuild (which only happens on
                              the first call when positions are static), the sort fires during
                              warmup. Use numWarmupCalls=0 to capture sort+NL rebuild cost.
*/
// =================================================================================================
struct BenchmarkScenario
{
    // ---- Particle geometry --------------------------------------------------
    uint              numParticles = 1000;
    ParticleShapeType shapeType    = ParticleShapeType::BOX;
    /** Half-extents / radius (X=Y=Z for sphere) used to build the Convex object. */
    Vector3<double> particleSize = {0.05, 0.05, 0.05};
    /** Stored in the CSV for bookkeeping; does not rescale particleSize here. */
    double aspectRatio = 1.0;

    // ---- Domain -------------------------------------------------------------
    Vector3<double> domainMin = {0.0, 0.0, 0.0};
    Vector3<double> domainMax = {1.0, 1.0, 1.0};

    // ---- Collision detection pipeline ---------------------------------------
    /** GJK sub-algorithm: GJK = Johnson, GJK_SV = SignedVolume. */
    NarrowPhaseType narrowPhaseType = NarrowPhaseType::GJK;
    /** Enable warm-start (simplex-carry-over) acceleration inside GJK. */
    bool gjkAcceleration = false;
    /** Bounding-volume pre-filter: OFF = none, OBB = oriented bounding box SAT.
        Active only on the relative-transform path (useRelativeTransformations=true). */
    BoundingVolumeType boundingVolumeType = BoundingVolumeType::OFF;
    /** true  = compute per-pair relative pos/quat, then run GJK in local frame
        false = feed world-frame pos/quat directly to GJK (no pre-transform step) */
    bool useRelativeTransformations = true;
    /** true  = pre-build ShapeData arrays (vtable-free GPU path)
        false = use Convex<T> virtual dispatch (CPU default path) */
    bool usePrebuiltShapes = false;

    // ---- Neighbor list ------------------------------------------------------
    NeighborListType neighborListType = NeighborListType::LINKEDCELL;
    /** LinkedCell variant (only matters when neighborListType == LINKEDCELL).
        - HOST       : classic CPU linked cells
        - SORTBASED  : GPU sort-based linked cells
        - ATOMIC     : GPU atomic linked cells
        - ATOMICFIXED: GPU atomic linked cells with pre-sized 2-D array */
    LinkedCellType linkedCellType = LinkedCellType::HOST;
    /** 0 = never sort particles; N = sort every N neighbor-list rebuilds. */
    uint sortFrequency = 0;
    /** Adaptive-skin update frequency (0 = no adaptive skin, always rebuild). */
    uint updateFrequency = 0;
    /** Initial estimate of pairs-per-particle for buffer pre-allocation. */
    uint initialPairsPerParticle = 16;

    // ---- Run control --------------------------------------------------------
    /** run() calls executed before the timer is reset (amortises NL rebuild overhead). */
    uint numWarmupCalls = 1;
    /** run() calls whose stage times are accumulated; CSV stores per-call averages. */
    uint numMeasureCalls = 10;
    uint randomSeed      = 42;
};

#endif
