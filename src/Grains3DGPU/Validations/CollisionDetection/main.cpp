#include "CollisionDetectionBenchmark.hh"
#include <cstring>
#include <iostream>
#include <vector>

template <typename T>
void runBenchmarkWithPrecision(const BenchmarkConfig& config,
                               const std::string&     csvFilename,
                               bool                   appendToCSV)
{
    CollisionDetectionBenchmark<T> benchmark(config, csvFilename, appendToCSV);
    benchmark.runBenchmark();
}

void printUsage(const char* progName)
{
    std::cout
        << "Usage: " << progName << " [options]\n"
        << "Options:\n"
        << "  --particles <N>           Number of particles (required)\n"
        << "  --domain <size>           Domain size (cube side length, required)\n"
        << "  --shape <type>            Shape type: sphere, box, superquadric (required)\n"
        << "  --size <x> <y> <z>        Particle size (3 values, required)\n"
        << "  --aspect <ratio>          Aspect ratio (default: 1.0)\n"
        << "  --precision <type>        Precision: single, double, both (default: both)\n"
        << "  --platform <type>         Platform: cpu, gpu, both (default: both)\n"
        << "  --trials <N>              Number of trials (default: 3)\n"
        << "  --seed <N>                Random seed (default: 42)\n"
        << "  --csv <filename>          CSV output file (default: data/collision_benchmark.csv)\n"
        << "  --append                  Append to CSV instead of overwriting\n"
        << "  --validate                Write contact info for validation\n"
        << "  --help                    Show this help message\n";
}

int main(int argc, char* argv[])
{
    // Parse command-line arguments
    uint              numParticles     = 0;
    double            domainSize       = 0.0;
    ParticleShapeType shapeType        = ParticleShapeType::SPHERE;
    Vector3<double>   particleSize     = Vector3<double>(0.0, 0.0, 0.0);
    double            aspectRatio      = 1.0;
    std::string       precisionStr     = "both";
    std::string       platformStr      = "both";
    uint              numTrials        = 3;
    uint              randomSeed       = 42;
    std::string       csvFilename      = "data/collision_benchmark.csv";
    bool              appendToCSV      = false;
    bool              validateContacts = false;

    bool hasParticles = false;
    bool hasDomain    = false;
    bool hasShape     = false;
    bool hasSize      = false;

    for(int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if(arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if(arg == "--particles" && i + 1 < argc)
        {
            numParticles = std::stoul(argv[++i]);
            hasParticles = true;
        }
        else if(arg == "--domain" && i + 1 < argc)
        {
            domainSize = std::stod(argv[++i]);
            hasDomain  = true;
        }
        else if(arg == "--shape" && i + 1 < argc)
        {
            std::string shape = argv[++i];
            if(shape == "sphere")
                shapeType = ParticleShapeType::SPHERE;
            else if(shape == "box")
                shapeType = ParticleShapeType::BOX;
            else if(shape == "superquadric")
                shapeType = ParticleShapeType::SUPERQUADRIC;
            else
            {
                std::cerr << "Error: Invalid shape type: " << shape << "\n";
                return 1;
            }
            hasShape = true;
        }
        else if(arg == "--size" && i + 3 < argc)
        {
            particleSize[0] = std::stod(argv[++i]);
            particleSize[1] = std::stod(argv[++i]);
            particleSize[2] = std::stod(argv[++i]);
            hasSize         = true;
        }
        else if(arg == "--aspect" && i + 1 < argc)
        {
            aspectRatio = std::stod(argv[++i]);
        }
        else if(arg == "--precision" && i + 1 < argc)
        {
            precisionStr = argv[++i];
        }
        else if(arg == "--platform" && i + 1 < argc)
        {
            platformStr = argv[++i];
        }
        else if(arg == "--trials" && i + 1 < argc)
        {
            numTrials = std::stoul(argv[++i]);
        }
        else if(arg == "--seed" && i + 1 < argc)
        {
            randomSeed = std::stoul(argv[++i]);
        }
        else if(arg == "--csv" && i + 1 < argc)
        {
            csvFilename = argv[++i];
        }
        else if(arg == "--append")
        {
            appendToCSV = true;
        }
        else if(arg == "--validate")
        {
            validateContacts = true;
        }
        else
        {
            std::cerr << "Error: Unknown argument or missing value: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Validate required arguments
    if(!hasParticles || !hasDomain || !hasShape || !hasSize)
    {
        std::cerr << "Error: Missing required arguments\n";
        printUsage(argv[0]);
        return 1;
    }

    // Configure domain
    Vector3<double> domainMin(0.0, 0.0, 0.0);
    Vector3<double> domainMax(domainSize, domainSize, domainSize);

    // Parse precision
    std::vector<PrecisionType> precisions;
    if(precisionStr == "single")
        precisions = {PrecisionType::SINGLE};
    else if(precisionStr == "double")
        precisions = {PrecisionType::DOUBLE};
    else if(precisionStr == "both")
        precisions = {PrecisionType::SINGLE, PrecisionType::DOUBLE};
    else
    {
        std::cerr << "Error: Invalid precision type: " << precisionStr << "\n";
        return 1;
    }

    // Parse platform
    PLATFORM platform;
    if(platformStr == "cpu")
        platform = PLATFORM::CPU;
    else if(platformStr == "gpu")
        platform = PLATFORM::GPU;
    else if(platformStr == "both")
        platform = PLATFORM::BOTH;
    else
    {
        std::cerr << "Error: Invalid platform type: " << platformStr << "\n";
        return 1;
    }

    // Run benchmarks for each precision
    int configCount = 0;

    // =========================================================================================
    // Iterate through precisions
    // =========================================================================================
    for(auto precision : precisions)
    {
        configCount++;

        // Create configuration
        BenchmarkConfig config;
        config.precision        = precision;
        config.platform         = platform;
        config.numParticles     = numParticles;
        config.shapeType        = shapeType;
        config.particleSize     = particleSize;
        config.aspectRatio      = aspectRatio;
        config.domainMin        = domainMin;
        config.domainMax        = domainMax;
        config.numTrials        = numTrials;
        config.randomSeed       = randomSeed;
        config.validateContacts = validateContacts;

        // Run benchmark
        if(config.precision == PrecisionType::SINGLE)
        {
            runBenchmarkWithPrecision<float>(config, csvFilename, appendToCSV);
        }
        else
        {
            runBenchmarkWithPrecision<double>(config, csvFilename, appendToCSV);
        }

        appendToCSV = true;  // After first run, always append
    }

    return 0;
}