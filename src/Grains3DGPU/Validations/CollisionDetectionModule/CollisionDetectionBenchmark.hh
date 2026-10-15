#ifndef _COLLISIONDETECTIONBENCHMARK_HH_
#define _COLLISIONDETECTIONBENCHMARK_HH_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "BenchmarkConfig.hh"
#include "Box.hh"
#include "CSVWriter.hh"
#include "CollisionDetectionModule.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "Insertion.hh"
#include "InsertionWindow.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "RigidBodyFactory.hh"
#include "Sphere.hh"
#include "Superquadric.hh"
#include "Torce.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief Reusable host-side particle state produced by prepareParticleData.
    Holds positions, orientations, and kinematics independent of LC type or platform. */
// =================================================================================================
template <typename T>
struct HostParticleData
{
    GrainsMemBuffer<Vector3<T>, MemType::HOST>    pos;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST> quat;
    GrainsMemBuffer<Kinematics<T>, MemType::HOST> kin;
};

// =================================================================================================
/** @brief Reusable device-side RigidBody pointer array.
    Created once per scenario and shared across all GPU runners with different LC types.
    ~DeviceParticleData calls freeDevice to clean up device heap allocations. */
// =================================================================================================
template <typename T>
struct DeviceParticleData
{
    GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE> rb;

    DeviceParticleData()                                     = default;
    DeviceParticleData(const DeviceParticleData&)            = delete;
    DeviceParticleData& operator=(const DeviceParticleData&) = delete;
    DeviceParticleData(DeviceParticleData&&)                 = default;
    ~DeviceParticleData()
    {
        RigidBodyFactory<T>::freeDevice(rb);
    }
};

// =================================================================================================
/** @brief Executes one BenchmarkScenario using CollisionDetectionModule.
    Templated on floating-point precision T and memory location M (HOST / DEVICE).

    ### Pipeline
      1. Build CollisionDetectionParameters from the scenario.
      2. Set GrainsParameters<T>::m_collisionDetection (static; read by CDM at runtime).
      3. Insert particles randomly into the domain.
      4. For M == DEVICE  : copy rigid-body pointer arrays to device via RigidBodyFactory.
         For M == HOST    : use the host buffers directly.
      5. Construct CollisionDetectionModule<T, M>.
      6. Warmup: call run() numWarmupCalls times (timer disabled).
      7. Measure: enable CDM timer, reset it, call run() numMeasureCalls times.
      8. Divide accumulated stage times by numMeasureCalls -> per-call averages.
      9. Write one CSV row.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, MemType M>
class BenchmarkRunner
{
public:
    // -------------------------------------------------------------------------
    /** @brief Constructor -- allocates particle data, inserts them, and builds the CDM.
        @param scenario  Full benchmark configuration
        @param csv       CSV writer to append the result row to */
    BenchmarkRunner(const BenchmarkScenario& scenario, CSVWriter& csv)
        : m_scenario(scenario)
        , m_csv(csv)
    {
        initCommon(scenario, csv);
        insertParticles(buildCDParams().linkedCellParameters);
        buildDeviceData();
    }

    // -------------------------------------------------------------------------
    /** @brief Constructor -- reuses pre-inserted host particle data, skipping insertion.
        Used for HOST (CPU) runners where no device data is needed.
        @param scenario  Full benchmark configuration
        @param csv       CSV writer to append the result row to
        @param pd        Pre-built host particle arrays (pos / quat / kin) */
    BenchmarkRunner(const BenchmarkScenario&   scenario,
                    CSVWriter&                 csv,
                    const HostParticleData<T>& pd)
        : m_scenario(scenario)
        , m_csv(csv)
        , m_ownsDeviceObjects(false)
    {
        initCommon(scenario, csv);
        const uint N = scenario.numParticles;
        for(uint i = 0; i < N; ++i)
        {
            m_h_pos[i]  = pd.pos[i];
            m_h_quat[i] = pd.quat[i];
            m_h_kin[i]  = pd.kin[i];
        }
        buildDeviceData();
    }

    // -------------------------------------------------------------------------
    /** @brief Constructor -- reuses pre-inserted host particle data AND the shared device RB
        pointer array, skipping both insertion and copyHostToDevice.
        The device RigidBody objects are owned by dpd and must outlive this runner.
        @param scenario  Full benchmark configuration
        @param csv       CSV writer to append the result row to
        @param pd        Pre-built host particle arrays (pos / quat / kin)
        @param dpd       Shared device RigidBody pointer array (not owned by this runner) */
    BenchmarkRunner(const BenchmarkScenario&   scenario,
                    CSVWriter&                 csv,
                    const HostParticleData<T>& pd,
                    DeviceParticleData<T>&     dpd)
        : m_scenario(scenario)
        , m_csv(csv)
        , m_ownsDeviceObjects(false)
    {
        initCommon(scenario, csv);
        const uint N = scenario.numParticles;
        for(uint i = 0; i < N; ++i)
        {
            m_h_pos[i]  = pd.pos[i];
            m_h_quat[i] = pd.quat[i];
            m_h_kin[i]  = pd.kin[i];
        }
        buildDeviceDataBorrowed(dpd);
    }

    // -------------------------------------------------------------------------
    /** @brief Build and return a DeviceParticleData for a given scenario.
        Performs copyHostToDevice once; the result is shared across all GPU runners
        that differ only in LC type, avoiding repeated device heap alloc/free cycles.
        @param scenario  Benchmark scenario */
    static DeviceParticleData<T> prepareDeviceData(const BenchmarkScenario& scenario)
    {
        cudaErrCheck(cudaGetDeviceProperties(&GrainsParameters<T>::m_GPU, 0));

        const T                    sx = static_cast<T>(scenario.particleSize[0]);
        const T                    sy = static_cast<T>(scenario.particleSize[1]);
        const T                    sz = static_cast<T>(scenario.particleSize[2]);
        std::unique_ptr<Convex<T>> shape;
        switch(scenario.shapeType)
        {
        case ParticleShapeType::SPHERE:
            shape = std::make_unique<Sphere<T>>(sx);
            break;
        case ParticleShapeType::BOX:
            shape = std::make_unique<Box<T>>(sx, sy, sz);
            break;
        case ParticleShapeType::SUPERQUADRIC:
            shape = std::make_unique<Superquadric<T>>(sx, sy, sz, T(4.0), T(4.0));
            break;
        }

        const uint                                    N = scenario.numParticles;
        GrainsMemBuffer<RigidBody<T>*, MemType::HOST> h_rb;
        h_rb.initialize(N);
        for(uint i = 0; i < N; ++i)
            h_rb[i] = new RigidBody<T>(shape->clone(), T(0), 0, 1);

        DeviceParticleData<T> dpd;
        RigidBodyFactory<T>::copyHostToDevice(h_rb, dpd.rb);

        for(uint i = 0; i < N; ++i)
            delete h_rb[i];

        return dpd;
    }
    static HostParticleData<T> prepareParticleData(const BenchmarkScenario& scenario)
    {
        const T                    sx = static_cast<T>(scenario.particleSize[0]);
        const T                    sy = static_cast<T>(scenario.particleSize[1]);
        const T                    sz = static_cast<T>(scenario.particleSize[2]);
        std::unique_ptr<Convex<T>> shape;
        switch(scenario.shapeType)
        {
        case ParticleShapeType::SPHERE:
            shape = std::make_unique<Sphere<T>>(sx);
            break;
        case ParticleShapeType::BOX:
            shape = std::make_unique<Box<T>>(sx, sy, sz);
            break;
        case ParticleShapeType::SUPERQUADRIC:
            shape = std::make_unique<Superquadric<T>>(sx, sy, sz, T(4.0), T(4.0));
            break;
        }

        const uint                                    N = scenario.numParticles;
        GrainsMemBuffer<RigidBody<T>*, MemType::HOST> h_rb;
        h_rb.initialize(N);
        for(uint i = 0; i < N; ++i)
            h_rb[i] = new RigidBody<T>(shape->clone(), T(0), 0, 1);

        HostParticleData<T> pd;
        pd.pos.initialize(N);
        pd.quat.initialize(N);
        pd.kin.initialize(N);
        for(uint i = 0; i < N; ++i)
        {
            pd.pos[i]  = {};
            pd.quat[i] = Quaternion<T>(T(1), T(0), T(0), T(0));
        }

        // Build minimal LinkedCellParameters just for insertion
        // (lc type does not affect particle placement)
        const T                 circR = shape->computeCircumscribedRadius();
        LinkedCellParameters<T> lcp;
        lcp.type                            = LinkedCellType::HOST;
        lcp.minCorner                       = Vector3<T>(static_cast<T>(scenario.domainMin[0]),
                                   static_cast<T>(scenario.domainMin[1]),
                                   static_cast<T>(scenario.domainMin[2]));
        lcp.maxCorner                       = Vector3<T>(static_cast<T>(scenario.domainMax[0]),
                                   static_cast<T>(scenario.domainMax[1]),
                                   static_cast<T>(scenario.domainMax[2]));
        lcp.minCellSize                     = T(2) * circR;
        lcp.cellSizeFactor                  = T(1);
        lcp.updateFrequency                 = scenario.updateFrequency;
        lcp.sortFrequency                   = scenario.sortFrequency;
        lcp.initialNumberOfPairsPerParticle = scenario.initialPairsPerParticle;

        Vector3<T>                      omin(T(0), T(0), T(0));
        Vector3<T>                      omax(T(2.0 * M_PI), T(2.0 * M_PI), T(2.0 * M_PI));
        std::vector<InsertionWindow<T>> posWin, oriWin;
        posWin.emplace_back(lcp.minCorner, lcp.maxCorner, scenario.randomSeed);
        oriWin.emplace_back(omin, omax, scenario.randomSeed + 1);

        auto ins = std::make_unique<Insertion<T>>();
        ins->setPositionInsertionInfo(posWin);
        ins->setOrientationInsertionInfo(oriWin);
        GrainsMemBuffer<uint, MemType::HOST>          bodyTag;
        GrainsMemBuffer<Vector3<T>, MemType::HOST>    localPos;
        GrainsMemBuffer<Quaternion<T>, MemType::HOST> localQuat;
        bodyTag.initialize(N);
        localPos.initialize(N);
        localQuat.initialize(N);
        for(uint i = 0; i < N; ++i)
        {
            bodyTag[i]   = 0u;
            localPos[i]  = Vector3<T>();
            localQuat[i] = Quaternion<T>(T(1), T(0), T(0), T(0));
        }
        ins->insert(&h_rb, pd.pos, pd.quat, pd.kin, lcp, 0u, N, bodyTag, localPos, localQuat);

        for(uint i = 0; i < N; ++i)
            delete h_rb[i];

        return pd;
    }

    // -------------------------------------------------------------------------
    /** @brief Destructor */
    ~BenchmarkRunner()
    {
        for(uint i = 0; i < m_scenario.numParticles; ++i)
            delete m_h_rb[i];
        if constexpr(M == MemType::DEVICE)
            if(m_ownsDeviceObjects)
                RigidBodyFactory<T>::freeDevice(m_rb);
    }

    // -------------------------------------------------------------------------
    /** @brief Execute warmup, then measure, then write one CSV row.
        @param runID  Row identifier written to the CSV */
    void run(uint runID)
    {
        using GP = GrainsParameters<T>;

        // Re-apply the full CD params each call (in case multiple runners share the process)
        GP::m_collisionDetection = buildCDParams();
        GP::m_cdmTimer.disable();  // timer OFF during warmup
        GP::m_simulationState = SimulationState<T>{};

        // ------------------------------------------------------------------
        // Warmup phase (timer OFF)
        // ------------------------------------------------------------------
        for(uint w = 0; w < m_scenario.numWarmupCalls; ++w)
            m_cdm->run(m_rb.getData(),
                       m_pos,
                       m_quat,
                       m_kin,
                       m_torce,
                       m_rbIds,
                       m_localPos,
                       m_localQuat,
                       m_masterSlot,
                       m_contactInfo,
                       m_pairList,
                       m_counts);

        // ------------------------------------------------------------------
        // Measurement phase (timer ON, accumulate over numMeasureCalls)
        // ------------------------------------------------------------------
        GP::m_cdmTimer.enable(M == MemType::DEVICE);  // timer ON for measurement
        GP::m_cdmTimer.reset();                       // zero all stage accumulators

        for(uint c = 0; c < m_scenario.numMeasureCalls; ++c)
            m_cdm->run(m_rb.getData(),
                       m_pos,
                       m_quat,
                       m_kin,
                       m_torce,
                       m_rbIds,
                       m_localPos,
                       m_localQuat,
                       m_masterSlot,
                       m_contactInfo,
                       m_pairList,
                       m_counts);

        if constexpr(M == MemType::DEVICE)
            cudaDeviceSynchronize();

        // ------------------------------------------------------------------
        // Compute per-call averages (seconds -> milliseconds)
        // ------------------------------------------------------------------
        const double inv = 1000.0 / static_cast<double>(m_scenario.numMeasureCalls);

        const double sortMs  = GP::m_cdmTimer[CDMStage::Sort] * inv;
        const double nlMs    = GP::m_cdmTimer[CDMStage::NeighborList] * inv;
        const double relMs   = GP::m_cdmTimer[CDMStage::RelativeTransform] * inv;
        const double bvMs    = GP::m_cdmTimer[CDMStage::BVFilter] * inv;
        const double gjkMs   = GP::m_cdmTimer[CDMStage::NarrowPhase] * inv;
        const double trnMs   = GP::m_cdmTimer[CDMStage::Transform] * inv;
        const double totalMs = sortMs + nlMs + relMs + bvMs + gjkMs + trnMs;

        // ------------------------------------------------------------------
        // Write CSV row
        // ------------------------------------------------------------------
        m_csv.writeRow(runID,
                       M == MemType::HOST ? "CPU" : "GPU",
                       sizeof(T) == sizeof(float) ? "single" : "double",
                       m_scenario.numParticles,
                       shapeStr(),
                       m_scenario.particleSize[0],
                       m_scenario.particleSize[1],
                       m_scenario.particleSize[2],
                       m_scenario.aspectRatio,
                       m_scenario.narrowPhaseType == NarrowPhaseType::GJK ? "johnson"
                                                                          : "signedvolume",
                       m_scenario.gjkAcceleration ? 1 : 0,
                       bvStr(),
                       m_scenario.useRelativeTransformations ? 1 : 0,
                       m_scenario.usePrebuiltShapes ? 1 : 0,
                       nlStr(),
                       lcStr(),
                       m_scenario.sortFrequency,
                       m_scenario.numWarmupCalls,
                       m_scenario.numMeasureCalls,
                       m_counts.numPairs,
                       sortMs,
                       nlMs,
                       relMs,
                       bvMs,
                       gjkMs,
                       trnMs,
                       totalMs);
    }

    /**
     * @brief Execute a Cartesian product over the provided parameter lists.
     *
     * Each combination updates `m_scenario` in-place and calls `run(runID++)`.
     * This avoids reconstructing the BenchmarkRunner (and the expensive particle
     * insertion / CDM construction) for every small configuration change.
     */
    void runCartesianProduct(uint&                                  runID,
                             const std::vector<NarrowPhaseType>&    narrowPhases,
                             const std::vector<bool>&               gjkAccels,
                             const std::vector<BoundingVolumeType>& bvs,
                             const std::vector<bool>&               relTransforms,
                             const std::vector<bool>&               prebuiltShapes,
                             const std::vector<NeighborListType>&   nls,
                             const std::vector<LinkedCellType>&     lcts,
                             const std::vector<uint>&               sortFreqs)
    {
        for(auto np : narrowPhases)
            for(auto accel : gjkAccels)
                for(auto bv : bvs)
                    for(auto rel : relTransforms)
                        for(auto prebuilt : prebuiltShapes)
                            for(auto nl : nls)
                                for(auto lc : lcts)
                                    for(auto sf : sortFreqs)
                                    {
                                        m_scenario.narrowPhaseType            = np;
                                        m_scenario.gjkAcceleration            = accel;
                                        m_scenario.boundingVolumeType         = bv;
                                        m_scenario.useRelativeTransformations = rel;
                                        m_scenario.usePrebuiltShapes          = prebuilt;
                                        m_scenario.neighborListType           = nl;
                                        m_scenario.linkedCellType             = lc;
                                        m_scenario.sortFrequency              = sf;

                                        run(runID++);
                                    }
    }

    // -------------------------------------------------------------------------
    /** @brief Write the CSV header row (call once before all runs). */
    static void writeHeader(CSVWriter& csv)
    {
        csv.writeHeader({"RunID",
                         "Platform",
                         "Precision",
                         "NParticles",
                         "Shape",
                         "SizeX",
                         "SizeY",
                         "SizeZ",
                         "AspectRatio",
                         "GJKVariant",
                         "GJKAccel",
                         "BoundingVolume",
                         "UseRelTransform",
                         "UsePrebuiltShapes",
                         "NeighborListType",
                         "LinkedCellType",
                         "SortFrequency",
                         "NumWarmupCalls",
                         "NumMeasureCalls",
                         "PairCount",
                         "SortTime_ms",
                         "NeighborListTime_ms",
                         "RelTransformTime_ms",
                         "BVFilterTime_ms",
                         "NarrowPhaseTime_ms",
                         "TransformTime_ms",
                         "TotalTime_ms"});
    }

private:
    // =========================================================================
    // Private helpers
    // =========================================================================

    /** @brief Shared init: GPU props, shape, CD params, host rb/pos/quat/kin alloc. */
    void initCommon(const BenchmarkScenario& scenario, CSVWriter&)
    {
        // 0. GPU device properties
        if constexpr(M == MemType::DEVICE)
            cudaErrCheck(cudaGetDeviceProperties(&GrainsParameters<T>::m_GPU, 0));

        // 1. Build convex shape
        const T sx = static_cast<T>(scenario.particleSize[0]);
        const T sy = static_cast<T>(scenario.particleSize[1]);
        const T sz = static_cast<T>(scenario.particleSize[2]);
        switch(scenario.shapeType)
        {
        case ParticleShapeType::SPHERE:
            m_shape = std::make_unique<Sphere<T>>(sx);
            break;
        case ParticleShapeType::BOX:
            m_shape = std::make_unique<Box<T>>(sx, sy, sz);
            break;
        case ParticleShapeType::SUPERQUADRIC:
            m_shape = std::make_unique<Superquadric<T>>(sx, sy, sz, T(4.0), T(4.0));
            break;
        }

        // 2. Collision detection parameters
        CollisionDetectionParameters<T> cdp       = buildCDParams();
        GrainsParameters<T>::m_collisionDetection = cdp;
        GrainsParameters<T>::m_simulationState    = SimulationState<T>{};

        // 3. Allocate HOST-side particle arrays
        const uint N = scenario.numParticles;
        m_h_rb.initialize(N);
        m_h_pos.initialize(N);
        m_h_quat.initialize(N);
        m_h_kin.initialize(N);
        m_h_bodyTags.initialize(N);
        m_h_localPos.initialize(N);
        m_h_localQuat.initialize(N);
        for(uint i = 0; i < N; ++i)
        {
            m_h_rb[i]        = new RigidBody<T>(m_shape->clone(), T(0), 0, 1);
            m_h_pos[i]       = {};
            m_h_quat[i]      = Quaternion<T>(T(1), T(0), T(0), T(0));
            m_h_bodyTags[i]  = 0u;
            m_h_localPos[i]  = Vector3<T>();
            m_h_localQuat[i] = Quaternion<T>(T(1), T(0), T(0), T(0));
        }
    }

    /** @brief Upload host data to M-typed buffers and construct the CDM (owned device RB). */
    void buildDeviceData()
    {
        const uint N   = m_scenario.numParticles;
        const auto cdp = buildCDParams();

        m_rb.initialize(N);
        m_pos.initialize(N);
        m_quat.initialize(N);
        m_kin.initialize(N);
        m_torce.initialize(N);
        m_rbIds.initialize(N);
        m_localPos.initialize(N);
        m_localQuat.initialize(N);
        m_masterSlot.initialize(N);

        if constexpr(M == MemType::HOST)
        {
            for(uint i = 0; i < N; ++i)
            {
                m_rb[i]         = m_h_rb[i];
                m_pos[i]        = m_h_pos[i];
                m_quat[i]       = m_h_quat[i];
                m_kin[i]        = m_h_kin[i];
                m_rbIds[i]      = 0u;
                m_localPos[i]   = Vector3<T>();
                m_localQuat[i]  = Quaternion<T>(T(1), T(0), T(0), T(0));
                m_masterSlot[i] = 0u;
            }
        }
        else
        {
            GrainsMemBuffer<uint, MemType::HOST>          h_rbIds;
            GrainsMemBuffer<Vector3<T>, MemType::HOST>    h_localPos;
            GrainsMemBuffer<Quaternion<T>, MemType::HOST> h_localQuat;
            GrainsMemBuffer<uint, MemType::HOST>          h_masterSlot;
            h_rbIds.initialize(N);
            h_localPos.initialize(N);
            h_localQuat.initialize(N);
            h_masterSlot.initialize(N);
            for(uint i = 0; i < N; ++i)
            {
                h_rbIds[i]      = 0u;
                h_localPos[i]   = Vector3<T>();
                h_localQuat[i]  = Quaternion<T>(T(1), T(0), T(0), T(0));
                h_masterSlot[i] = 0u;
            }
            m_rbIds.copyFrom(h_rbIds);
            m_localPos.copyFrom(h_localPos);
            m_localQuat.copyFrom(h_localQuat);
            m_masterSlot.copyFrom(h_masterSlot);
            m_pos.copyFrom(m_h_pos);
            m_quat.copyFrom(m_h_quat);
            m_kin.copyFrom(m_h_kin);
            RigidBodyFactory<T>::copyHostToDevice(m_h_rb, m_rb);
        }

        const uint estPairs = m_scenario.initialPairsPerParticle * N;
        m_contactInfo.initialize(estPairs);
        m_pairList.initialize(estPairs);
        m_counts.numParticles = N;

        m_cdm = std::make_unique<CollisionDetectionModule<T, M>>(&m_rb,
                                                                 m_pos,
                                                                 m_quat,
                                                                 m_rbIds,
                                                                 m_h_bodyTags.getData(),
                                                                 cdp,
                                                                 0u,
                                                                 N);
    }

    /** @brief Variant of buildDeviceData that borrows the device RB pointer array from dpd
        instead of calling copyHostToDevice.  The actual RigidBody objects on the device heap
        remain owned by dpd. */
    void buildDeviceDataBorrowed(DeviceParticleData<T>& dpd)
    {
        static_assert(M == MemType::DEVICE,
                      "buildDeviceDataBorrowed is only valid for DEVICE runners");

        const uint N   = m_scenario.numParticles;
        const auto cdp = buildCDParams();

        // Allocate a new pointer array and copy device-side pointers from the shared buffer.
        // The pointed-to RigidBody objects are NOT duplicated -- only the pointer values are
        // copied.
        m_rb.initialize(N);
        cudaErrCheck(cudaMemcpy(m_rb.getData(),
                                dpd.rb.getData(),
                                N * sizeof(RigidBody<T>*),
                                cudaMemcpyDeviceToDevice));

        m_pos.initialize(N);
        m_quat.initialize(N);
        m_kin.initialize(N);
        m_torce.initialize(N);
        m_rbIds.initialize(N);
        m_localPos.initialize(N);
        m_localQuat.initialize(N);
        m_masterSlot.initialize(N);

        GrainsMemBuffer<uint, MemType::HOST>          h_rbIds;
        GrainsMemBuffer<Vector3<T>, MemType::HOST>    h_localPos;
        GrainsMemBuffer<Quaternion<T>, MemType::HOST> h_localQuat;
        GrainsMemBuffer<uint, MemType::HOST>          h_masterSlot;
        h_rbIds.initialize(N);
        h_localPos.initialize(N);
        h_localQuat.initialize(N);
        h_masterSlot.initialize(N);
        for(uint i = 0; i < N; ++i)
        {
            h_rbIds[i]      = 0u;
            h_localPos[i]   = Vector3<T>();
            h_localQuat[i]  = Quaternion<T>(T(1), T(0), T(0), T(0));
            h_masterSlot[i] = 0u;
        }
        m_rbIds.copyFrom(h_rbIds);
        m_localPos.copyFrom(h_localPos);
        m_localQuat.copyFrom(h_localQuat);
        m_masterSlot.copyFrom(h_masterSlot);
        m_pos.copyFrom(m_h_pos);
        m_quat.copyFrom(m_h_quat);
        m_kin.copyFrom(m_h_kin);

        const uint estPairs = m_scenario.initialPairsPerParticle * N;
        m_contactInfo.initialize(estPairs);
        m_pairList.initialize(estPairs);
        m_counts.numParticles = N;

        m_cdm = std::make_unique<CollisionDetectionModule<T, M>>(&m_rb,
                                                                 m_pos,
                                                                 m_quat,
                                                                 m_rbIds,
                                                                 m_h_bodyTags.getData(),
                                                                 cdp,
                                                                 0u,
                                                                 N);
    }

    /** @brief Build a CollisionDetectionParameters from the current scenario. */
    CollisionDetectionParameters<T> buildCDParams() const
    {
        CollisionDetectionParameters<T> cdp;
        cdp.neighborListType           = m_scenario.neighborListType;
        cdp.boundingVolumeType         = m_scenario.boundingVolumeType;
        cdp.narrowPhaseType            = m_scenario.narrowPhaseType;
        cdp.gjkAcceleration            = m_scenario.gjkAcceleration;
        cdp.useRelativeTransformations = m_scenario.useRelativeTransformations;
        cdp.usePrebuiltShapes          = m_scenario.usePrebuiltShapes;

        LinkedCellParameters<T>& lcp = cdp.linkedCellParameters;
        lcp.type                     = m_scenario.linkedCellType;
        lcp.minCorner                = Vector3<T>(static_cast<T>(m_scenario.domainMin[0]),
                                   static_cast<T>(m_scenario.domainMin[1]),
                                   static_cast<T>(m_scenario.domainMin[2]));
        lcp.maxCorner                = Vector3<T>(static_cast<T>(m_scenario.domainMax[0]),
                                   static_cast<T>(m_scenario.domainMax[1]),
                                   static_cast<T>(m_scenario.domainMax[2]));
        // compute circumscribed radius from the constructed shape if available;
        // otherwise fall back to a conservative estimate from particleSize.
        T circR;
        if(m_shape)
            circR = m_shape->computeCircumscribedRadius();
        else
        {
            const T sx = static_cast<T>(m_scenario.particleSize[0]);
            const T sy = static_cast<T>(m_scenario.particleSize[1]);
            const T sz = static_cast<T>(m_scenario.particleSize[2]);
            circR      = std::sqrt(sx * sx + sy * sy + sz * sz);
        }
        lcp.minCellSize                     = T(2) * circR;
        lcp.cellSizeFactor                  = T(1);
        lcp.updateFrequency                 = m_scenario.updateFrequency;
        lcp.sortFrequency                   = m_scenario.sortFrequency;
        lcp.initialNumberOfPairsPerParticle = m_scenario.initialPairsPerParticle;
        return cdp;
    }

    /** @brief Place particles randomly using the Insertion class. */
    void insertParticles(const LinkedCellParameters<T>& lcp)
    {
        const uint N = m_scenario.numParticles;
        Vector3<T> pmin(static_cast<T>(m_scenario.domainMin[0]),
                        static_cast<T>(m_scenario.domainMin[1]),
                        static_cast<T>(m_scenario.domainMin[2]));
        Vector3<T> pmax(static_cast<T>(m_scenario.domainMax[0]),
                        static_cast<T>(m_scenario.domainMax[1]),
                        static_cast<T>(m_scenario.domainMax[2]));
        Vector3<T> omin(T(0), T(0), T(0));
        Vector3<T> omax(T(2.0 * M_PI), T(2.0 * M_PI), T(2.0 * M_PI));

        std::vector<InsertionWindow<T>> posWin;
        posWin.emplace_back(pmin, pmax, m_scenario.randomSeed);
        std::vector<InsertionWindow<T>> oriWin;
        oriWin.emplace_back(omin, omax, m_scenario.randomSeed + 1);

        auto ins = std::make_unique<Insertion<T>>();
        ins->setPositionInsertionInfo(posWin);
        ins->setOrientationInsertionInfo(oriWin);
        ins->insert(&m_h_rb,
                    m_h_pos,
                    m_h_quat,
                    m_h_kin,
                    lcp,
                    0u,
                    N,
                    m_h_bodyTags,
                    m_h_localPos,
                    m_h_localQuat);
    }

    // -- Label helpers for CSV strings ----------------------------------------
    const char* shapeStr() const
    {
        switch(m_scenario.shapeType)
        {
        case ParticleShapeType::SPHERE:
            return "sphere";
        case ParticleShapeType::BOX:
            return "box";
        case ParticleShapeType::SUPERQUADRIC:
            return "superquadric";
        }
        return "unknown";
    }
    const char* bvStr() const
    {
        switch(m_scenario.boundingVolumeType)
        {
        case BoundingVolumeType::OFF:
            return "off";
        case BoundingVolumeType::OBB:
            return "obb";
        case BoundingVolumeType::OBC:
            return "obc";
        }
        return "off";
    }
    const char* nlStr() const
    {
        return m_scenario.neighborListType == NeighborListType::LINKEDCELL ? "linkedcell" : "nsq";
    }
    const char* lcStr() const
    {
        switch(m_scenario.linkedCellType)
        {
        case LinkedCellType::HOST:
            return "host";
        case LinkedCellType::SORTBASED:
            return "sortbased";
        case LinkedCellType::ATOMIC:
            return "atomic";
        case LinkedCellType::ATOMICFIXED:
            return "atomicfixed";
        }
        return "host";
    }

    // =========================================================================
    // Data members
    // =========================================================================
    BenchmarkScenario m_scenario;
    CSVWriter&        m_csv;
    /** True when this runner created device RigidBody objects (must call freeDevice).
        False when device objects are shared via DeviceParticleData (borrowed). */
    bool m_ownsDeviceObjects = true;

    std::unique_ptr<Convex<T>> m_shape;  ///< Shape template (owned here)

    // HOST-side particle arrays (used for insertion + as a stable owner of RigidBody*)
    GrainsMemBuffer<RigidBody<T>*, MemType::HOST> m_h_rb;
    GrainsMemBuffer<Vector3<T>, MemType::HOST>    m_h_pos;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST> m_h_quat;
    GrainsMemBuffer<Kinematics<T>, MemType::HOST> m_h_kin;
    GrainsMemBuffer<uint, MemType::HOST>          m_h_bodyTags;   ///< shape-ID 0 for all particles
    GrainsMemBuffer<Vector3<T>, MemType::HOST>    m_h_localPos;   ///< local pos offsets (zero)
    GrainsMemBuffer<Quaternion<T>, MemType::HOST> m_h_localQuat;  ///< local quat offsets (identity)

    // M-typed arrays passed to CDM at construction and at every run() call
    GrainsMemBuffer<RigidBody<T>*, M> m_rb;
    GrainsMemBuffer<Vector3<T>, M>    m_pos;
    GrainsMemBuffer<Quaternion<T>, M> m_quat;
    GrainsMemBuffer<Kinematics<T>, M> m_kin;
    GrainsMemBuffer<Torce<T>, M>      m_torce;
    GrainsMemBuffer<uint, M>          m_rbIds;       ///< body tags (sequential, no composites)
    GrainsMemBuffer<Vector3<T>, M>    m_localPos;    ///< local pos offsets (all zero)
    GrainsMemBuffer<Quaternion<T>, M> m_localQuat;   ///< local quat offsets (all identity)
    GrainsMemBuffer<uint, M>          m_masterSlot;  ///< master slot lookup (no composites)

    // CDM-owned pair buffers resized by CDM internally; pre-allocated here
    GrainsMemBuffer<ContactInfo<T>, M> m_contactInfo;
    GrainsMemBuffer<uint2, M>          m_pairList;
    ComponentCounts                    m_counts;  ///< numParticles set in buildDeviceData

    std::unique_ptr<CollisionDetectionModule<T, M>> m_cdm;
};

#endif
