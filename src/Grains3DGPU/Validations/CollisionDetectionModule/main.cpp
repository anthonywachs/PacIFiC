#include "BenchmarkConfig.hh"
#include "CSVWriter.hh"
#include "CollisionDetectionBenchmark.hh"

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// =================================================================================================
// Forward declarations
// =================================================================================================

/** @brief Dispatch to BenchmarkRunner<T, M> for the right precision x platform combo.
    Runs the full Cartesian product over all CD-pipeline axes.
    runID is incremented in-place. */
template <typename T>
static void dispatchPlatform(const BenchmarkScenario&               s,
                             bool                                   runCPU,
                             bool                                   runGPU,
                             CSVWriter&                             csv,
                             uint&                                  runID,
                             const std::vector<NarrowPhaseType>&    gjkTypes,
                             const std::vector<bool>&               gjkAccels,
                             const std::vector<BoundingVolumeType>& bvTypes,
                             const std::vector<bool>&               relTransforms,
                             const std::vector<bool>&               prebuiltShapes,
                             const std::vector<LinkedCellType>&     cpuLcTypes,
                             const std::vector<LinkedCellType>&     gpuLcTypes,
                             const std::vector<uint>&               sortFreqs)
{
    // Insert particles once (independent of LC type / platform / precision).
    const BenchmarkScenario scBase = [&] {
        BenchmarkScenario sc = s;
        sc.linkedCellType    = cpuLcTypes[0];
        return sc;
    }();
    const HostParticleData<T> pd = BenchmarkRunner<T, MemType::HOST>::prepareParticleData(scBase);

    // NSQ stores all N*(N-1)/2 pairs; cap to avoid OOM and excessive narrow-phase work.
    // Only run NSQ for N <= 2048.
    static constexpr uint NSQ_MAX_N = 2048;
    const bool            canRunNSQ = (s.numParticles <= NSQ_MAX_N);

    if(runCPU)
    {
        // --- NSQ pass (once; LC type irrelevant) ---------------------------------
        // IMPORTANT: neighborListType must be set on the scenario BEFORE
        // constructing BenchmarkRunner so the CDM is built with the correct NL.
        if(canRunNSQ)
        {
            BenchmarkScenario sc = s;
            sc.neighborListType  = NeighborListType::NSQ;
            sc.linkedCellType    = cpuLcTypes[0];  // placeholder, unused by NSQ
            BenchmarkRunner<T, MemType::HOST>(sc, csv, pd)
                .runCartesianProduct(runID,
                                     gjkTypes,
                                     gjkAccels,
                                     bvTypes,
                                     relTransforms,
                                     prebuiltShapes,
                                     {NeighborListType::NSQ},
                                     {cpuLcTypes[0]},
                                     sortFreqs);
        }

        // --- LINKEDCELL pass (one runner per CPU LC type) ------------------------
        for(auto lct : cpuLcTypes)
        {
            BenchmarkScenario sc = s;
            sc.neighborListType  = NeighborListType::LINKEDCELL;
            sc.linkedCellType    = lct;
            BenchmarkRunner<T, MemType::HOST>(sc, csv, pd)
                .runCartesianProduct(runID,
                                     gjkTypes,
                                     gjkAccels,
                                     bvTypes,
                                     relTransforms,
                                     prebuiltShapes,
                                     {NeighborListType::LINKEDCELL},
                                     {lct},
                                     sortFreqs);
        }
    }
    if(runGPU)
    {
        // Build device RigidBody objects ONCE and share them across all LC-type runners.
        // This avoids repeated device heap alloc/free cycles which fragment the device heap.
        DeviceParticleData<T> dpd = BenchmarkRunner<T, MemType::HOST>::prepareDeviceData(scBase);

        // --- NSQ pass (once; LC type irrelevant) ---------------------------------
        if(canRunNSQ)
        {
            BenchmarkScenario sc = s;
            sc.neighborListType  = NeighborListType::NSQ;
            sc.linkedCellType    = gpuLcTypes[0];  // placeholder, unused by NSQ
            BenchmarkRunner<T, MemType::DEVICE>(sc, csv, pd, dpd)
                .runCartesianProduct(runID,
                                     gjkTypes,
                                     gjkAccels,
                                     bvTypes,
                                     relTransforms,
                                     prebuiltShapes,
                                     {NeighborListType::NSQ},
                                     {gpuLcTypes[0]},
                                     sortFreqs);
        }

        // --- LINKEDCELL pass (one runner per GPU LC type) ------------------------
        for(auto lct : gpuLcTypes)
        {
            BenchmarkScenario sc = s;
            sc.neighborListType  = NeighborListType::LINKEDCELL;
            sc.linkedCellType    = lct;
            BenchmarkRunner<T, MemType::DEVICE>(sc, csv, pd, dpd)
                .runCartesianProduct(runID,
                                     gjkTypes,
                                     gjkAccels,
                                     bvTypes,
                                     relTransforms,
                                     prebuiltShapes,
                                     {NeighborListType::LINKEDCELL},
                                     {lct},
                                     sortFreqs);
        }
        // dpd goes out of scope here -> frees device RigidBody objects once
    }
}

// =================================================================================================
// Usage
// =================================================================================================
static void printUsage(const char* prog)
{
    std::cout
        << "\nUsage: " << prog << " [options]\n\n"
        << "  All collision-detection parameter combinations are run automatically.\n\n"
        << "Required:\n"
        << "  --particles N           Number of particles\n"
        << "  --shape sphere|box|superquadric\n"
        << "  --size X Y Z            Half-extents (sphere: X=radius; box/SQ: full half-extents)\n"
        << "\nDomain:\n"
        << "  --domain S              Cubic domain [0,S]^3  (default: 1.0)\n"
        << "  --domain-min X Y Z      Domain lower corner\n"
        << "  --domain-max X Y Z      Domain upper corner\n"
        << "\nNeighbor list:\n"
        << "  --update-freq N              Adaptive skin update freq, 0=disabled (default: 0)\n"
        << "  --pairs-per-particle N       Initial pair buffer size hint (default: 16)\n"
        << "\nPrecision / platform:\n"
        << "  --precision single|double|both  (default: both)\n"
        << "  --platform cpu|gpu|both         (default: both)\n"
        << "  --prebuilt on|off|both          PrebuiltShapes vtable-free path (default: both)\n"
        << "\nRun control:\n"
        << "  --warmup N               Warmup calls before timing (default: 1)\n"
        << "  --measure N              Timed calls per scenario (default: 10)\n"
        << "  --seed N                 Random seed (default: 42)\n"
        << "  --aspect R               Aspect ratio label stored in CSV (default: 1.0)\n"
        << "\nOutput:\n"
        << "  --csv FILE               Output CSV path (default: data/benchmark.csv)\n"
        << "  --append                 Append to existing CSV instead of overwriting\n"
        << "  --help\n\n";
}

// =================================================================================================
// main
// =================================================================================================
int main(int argc, char* argv[])
{
    // ---- defaults -----------------------------------------------------------
    BenchmarkScenario s;
    s.numParticles            = 0;  // 0 = not set
    s.shapeType               = ParticleShapeType::BOX;
    s.particleSize            = {0.05, 0.05, 0.05};
    s.aspectRatio             = 1.0;
    s.domainMin               = {0.0, 0.0, 0.0};
    s.domainMax               = {1.0, 1.0, 1.0};
    s.updateFrequency         = 0;
    s.initialPairsPerParticle = 16;
    s.numWarmupCalls          = 1;
    s.numMeasureCalls         = 10;
    s.randomSeed              = 42;

    std::string precStr      = "both";
    std::string platStr      = "both";
    std::string prebuiltStr  = "both";
    std::string csvFile      = "data/benchmark.csv";
    bool        appendCSV    = false;
    bool        hasParticles = false;
    bool        hasShape     = false;
    bool        hasSize      = false;

    // ---- argument parsing ---------------------------------------------------
    for(int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

        if(a == "--help" || a == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if(a == "--particles" && i + 1 < argc)
        {
            s.numParticles = std::stoul(argv[++i]);
            hasParticles   = true;
        }
        else if(a == "--shape" && i + 1 < argc)
        {
            std::string sh = argv[++i];
            if(sh == "sphere")
                s.shapeType = ParticleShapeType::SPHERE;
            else if(sh == "box")
                s.shapeType = ParticleShapeType::BOX;
            else if(sh == "superquadric")
                s.shapeType = ParticleShapeType::SUPERQUADRIC;
            else
            {
                std::cerr << "Unknown shape: " << sh << "\n";
                return 1;
            }
            hasShape = true;
        }
        else if(a == "--size" && i + 3 < argc)
        {
            s.particleSize[0] = std::stod(argv[++i]);
            s.particleSize[1] = std::stod(argv[++i]);
            s.particleSize[2] = std::stod(argv[++i]);
            hasSize           = true;
        }
        else if(a == "--aspect" && i + 1 < argc)
        {
            s.aspectRatio = std::stod(argv[++i]);
        }
        else if(a == "--domain" && i + 1 < argc)
        {
            double d    = std::stod(argv[++i]);
            s.domainMin = {0.0, 0.0, 0.0};
            s.domainMax = {d, d, d};
        }
        else if(a == "--domain-min" && i + 3 < argc)
        {
            s.domainMin[0] = std::stod(argv[++i]);
            s.domainMin[1] = std::stod(argv[++i]);
            s.domainMin[2] = std::stod(argv[++i]);
        }
        else if(a == "--domain-max" && i + 3 < argc)
        {
            s.domainMax[0] = std::stod(argv[++i]);
            s.domainMax[1] = std::stod(argv[++i]);
            s.domainMax[2] = std::stod(argv[++i]);
        }
        else if(a == "--update-freq" && i + 1 < argc)
        {
            s.updateFrequency = std::stoul(argv[++i]);
        }
        else if(a == "--pairs-per-particle" && i + 1 < argc)
        {
            s.initialPairsPerParticle = std::stoul(argv[++i]);
        }
        else if(a == "--warmup" && i + 1 < argc)
        {
            s.numWarmupCalls = std::stoul(argv[++i]);
        }
        else if(a == "--measure" && i + 1 < argc)
        {
            s.numMeasureCalls = std::stoul(argv[++i]);
        }
        else if(a == "--seed" && i + 1 < argc)
        {
            s.randomSeed = std::stoul(argv[++i]);
        }
        else if(a == "--precision" && i + 1 < argc)
        {
            precStr = argv[++i];
        }
        else if(a == "--platform" && i + 1 < argc)
        {
            platStr = argv[++i];
        }
        else if(a == "--prebuilt" && i + 1 < argc)
        {
            prebuiltStr = argv[++i];
        }
        else if(a == "--csv" && i + 1 < argc)
        {
            csvFile = argv[++i];
        }
        else if(a == "--append")
        {
            appendCSV = true;
        }
        else
        {
            std::cerr << "Unknown argument: " << a << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // ---- validate -----------------------------------------------------------
    if(!hasParticles)
    {
        std::cerr << "Error: --particles is required\n";
        printUsage(argv[0]);
        return 1;
    }
    if(!hasShape)
    {
        std::cerr << "Error: --shape is required\n";
        return 1;
    }
    if(!hasSize)
    {
        std::cerr << "Error: --size is required\n";
        return 1;
    }

    // ---- all collision-detection combinations (hard-coded) ------------------
    const std::vector<NarrowPhaseType> gjkTypes  = {NarrowPhaseType::GJK, NarrowPhaseType::GJK_SV};
    const std::vector<bool>            gjkAccels = {false, true};
    const std::vector<BoundingVolumeType> bvTypes
        = {BoundingVolumeType::OFF, BoundingVolumeType::OBB, BoundingVolumeType::OBC};
    const std::vector<bool> relTransforms        = {false, true};
    const std::vector<bool> prebuiltShapes       = (prebuiltStr == "on") ? std::vector<bool>{true}
                                                   : (prebuiltStr == "off")
                                                       ? std::vector<bool>{false}
                                                       : std::vector<bool>{false, true};
    const std::vector<LinkedCellType> cpuLcTypes = {LinkedCellType::HOST};
    const std::vector<LinkedCellType> gpuLcTypes
        = {LinkedCellType::SORTBASED, LinkedCellType::ATOMIC, LinkedCellType::ATOMICFIXED};
    const std::vector<uint> sortFreqs = {0};

    // ---- precision / platform lists -----------------------------------------
    bool runSingle = (precStr == "single" || precStr == "both");
    bool runDouble = (precStr == "double" || precStr == "both");
    bool runCPU    = (platStr == "cpu" || platStr == "both");
    bool runGPU    = (platStr == "gpu" || platStr == "both");

    if(!runSingle && !runDouble)
    {
        std::cerr << "Unknown precision: " << precStr << "\n";
        return 1;
    }
    if(!runCPU && !runGPU)
    {
        std::cerr << "Unknown platform: " << platStr << "\n";
        return 1;
    }

    // ---- CSV setup ----------------------------------------------------------
    CSVWriter csv(csvFile, appendCSV);
    if(!appendCSV)
        BenchmarkRunner<float, MemType::HOST>::writeHeader(csv);

    uint runID = 0;

    // ---- dispatch -----------------------------------------------------------
    if(runSingle)
        dispatchPlatform<float>(s,
                                runCPU,
                                runGPU,
                                csv,
                                runID,
                                gjkTypes,
                                gjkAccels,
                                bvTypes,
                                relTransforms,
                                prebuiltShapes,
                                cpuLcTypes,
                                gpuLcTypes,
                                sortFreqs);
    if(runDouble)
        dispatchPlatform<double>(s,
                                 runCPU,
                                 runGPU,
                                 csv,
                                 runID,
                                 gjkTypes,
                                 gjkAccels,
                                 bvTypes,
                                 relTransforms,
                                 prebuiltShapes,
                                 cpuLcTypes,
                                 gpuLcTypes,
                                 sortFreqs);

    return 0;
}
