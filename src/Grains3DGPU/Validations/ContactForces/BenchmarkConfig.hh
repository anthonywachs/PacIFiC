#ifndef _BENCHMARKCONFIG_HH_
#define _BENCHMARKCONFIG_HH_

enum class PLATFORM
{
    CPU  = 0,
    GPU  = 1,
    BOTH = 2
};

// =================================================================================================
/** @brief Configuration for Contact Table benchmarks */
struct BenchmarkConfig
{
    // Platform parameters
    PLATFORM platform;

    // Parameters
    uint   capacity;
    double loadFactor;
    uint   numLookups;

    // Test parameters
    uint numTrials;
    uint randomSeed;
    bool validateContacts;

    BenchmarkConfig()
        : platform(PLATFORM::BOTH)
        , capacity(2000)
        , loadFactor(0.5)
        , numLookups(10000)
        , numTrials(5)
        , randomSeed(42)
        , validateContacts(false){};
};

#endif
