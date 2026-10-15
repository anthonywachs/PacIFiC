#ifndef _CONTACTTABLEBENCHMARK_HH_
#define _CONTACTTABLEBENCHMARK_HH_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <cuda_runtime.h>

#include "BenchmarkConfig.hh"
#include "CSVWriter.hh"
#include "ContactTable.hh"
#include "GrainsUtils.hh"
#include "StepTimer.hh"

// =================================================================================================
/** @brief Contact Hash Table Performance Benchmark

    This class benchmarks ContactHashTable performance for insert and lookup operations.
    Measures host serial performance and GPU kernel performance with proper pre-population.
    Results are exported to CSV for post-processing.

    @author A.Yazdani - 2026 - Contact Table Performance Validation */
// =================================================================================================

__global__ static void
    bench_findOrInsert_kernel(ContactHashTableView view, uint2* keys, uint* out, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= n)
        return;
    view.findOrInsert(keys[idx], out[idx]);
}

class ContactTableBenchmark
{
private:
    BenchmarkConfig            m_config;
    std::unique_ptr<CSVWriter> m_csvWriter;

    // Result structure for benchmark output
    struct Result
    {
        int    N;
        double host_ms;
        double kernel_ms;
        double pipeline_ms;
    };

public:
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with configuration */
    ContactTableBenchmark(const BenchmarkConfig& config,
                          const std::string&     csvFilename,
                          bool                   appendToCSV = false)
        : m_config(config)
    {
        m_csvWriter = std::make_unique<CSVWriter>(csvFilename, appendToCSV);

        // Initialize CSV header only if not appending
        if(!appendToCSV)
        {
            std::vector<std::string> columns = {"numTrials",
                                                "capacity",
                                                "loadFactor",
                                                "numLookups",
                                                "host_ms",
                                                "kernel_ms",
                                                "pipeline_ms"};
            m_csvWriter->writeHeader(columns);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~ContactTableBenchmark() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Run benchmark with current configuration */
    void runBenchmark()
    {
        Result result = runBenchmarkInternal(m_config);
        appendCSV(result);
    }

private:
    // ---------------------------------------------------------------------------------------------
    /** @brief Generate random key set
        @param N Number of keys to generate
        @param seed Random seed (0 = deterministic)
        @return Vector of uint2 keys */
    std::vector<uint2> makeKeys(int N, unsigned seed = 0)
    {
        std::vector<uint2> keys(N);
        for(int i = 0; i < N; i++)
        {
            keys[i].x = i;
            keys[i].y = 0;
        }
        if(seed != 0)
        {
            std::mt19937 rng(seed);
            std::shuffle(keys.begin(), keys.end(), rng);
        }
        return keys;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Run host benchmark once
        @param N Number of lookup operations
        @param keys Lookup keys
        @param N_preload Number of keys to pre-populate
        @param preload_keys Keys for pre-population
        @param capacity Table capacity
        @return Host execution time in milliseconds */
    double runHostOnce(int                       N,
                       const std::vector<uint2>& keys,
                       int                       N_preload,
                       const std::vector<uint2>& preload_keys,
                       int                       capacity)
    {
        ContactHashTable<MemType::HOST> table(capacity);

        // Pre-populate table
        for(int i = 0; i < N_preload; i++)
        {
            uint idx;
            table.findOrInsert(preload_keys[i], idx);
        }

        // Time lookup operations
        StepTimer timer;
        timer.start();
        for(int i = 0; i < N; i++)
        {
            uint idx;
            table.findOrInsert(keys[i], idx);
        }
        timer.stop();
        return timer.getElapsedMilliseconds();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Run device benchmark once
        @param N Number of lookup operations
        @param keys Lookup keys
        @param N_preload Number of keys to pre-populate
        @param preload_keys Keys for pre-population
        @param capacity Table capacity
        @return Pair of (kernel_ms, pipeline_ms) */
    std::pair<float, float> runDeviceOnce(int                       N,
                                          const std::vector<uint2>& keys,
                                          int                       N_preload,
                                          const std::vector<uint2>& preload_keys,
                                          int                       capacity)
    {
        ContactHashTable<MemType::DEVICE> devTable(capacity);
        ContactHashTableView              view = devTable.getView();

        uint2* d_keys         = nullptr;
        uint*  d_out          = nullptr;
        uint2* d_preload_keys = nullptr;
        uint*  d_preload_out  = nullptr;

        cudaMalloc(&d_keys, N * sizeof(uint2));
        cudaMalloc(&d_out, N * sizeof(uint));
        cudaMalloc(&d_preload_keys, N_preload * sizeof(uint2));
        cudaMalloc(&d_preload_out, N_preload * sizeof(uint));

        cudaMemcpy(d_keys, keys.data(), N * sizeof(uint2), cudaMemcpyHostToDevice);
        cudaMemcpy(d_preload_keys,
                   preload_keys.data(),
                   N_preload * sizeof(uint2),
                   cudaMemcpyHostToDevice);

        // Pre-populate the table with N_preload entries
        {
            int devId = 0;
            cudaGetDevice(&devId);
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, devId);
            uint numBlocksU = 0, numThreadsU = 0;
            computeOptimalThreadsAndBlocks(static_cast<uint>(N_preload),
                                           prop,
                                           numBlocksU,
                                           numThreadsU);
            int blockSize_preload = static_cast<int>(numThreadsU);
            int numBlocks_preload = static_cast<int>(numBlocksU);
            bench_findOrInsert_kernel<<<numBlocks_preload, blockSize_preload>>>(view,
                                                                                d_preload_keys,
                                                                                d_preload_out,
                                                                                N_preload);
            cudaDeviceSynchronize();
        }

        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        // Compute optimal threads and blocks for current device and N
        int devId = 0;
        cudaGetDevice(&devId);
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, devId);
        uint numBlocksU = 0, numThreadsU = 0;
        computeOptimalThreadsAndBlocks(static_cast<uint>(N), prop, numBlocksU, numThreadsU);
        int blockSize = static_cast<int>(numThreadsU);
        int numBlocks = static_cast<int>(numBlocksU);

        // Warmup kernel to eliminate first-run overhead
        bench_findOrInsert_kernel<<<numBlocks, blockSize>>>(view, d_keys, d_out, N);
        cudaDeviceSynchronize();

        // Kernel timing using CUDA events
        cudaEventRecord(start);
        bench_findOrInsert_kernel<<<numBlocks, blockSize>>>(view, d_keys, d_out, N);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float kernel_ms = 0.0f;
        cudaEventElapsedTime(&kernel_ms, start, stop);

        // Full pipeline timing (H2D + kernel + D2H) measured with CUDA events
        std::vector<uint> h_out(N);
        float             pipeline_ms = 0.0f;
        cudaEventRecord(start);
        cudaMemcpy(d_keys, keys.data(), N * sizeof(uint2), cudaMemcpyHostToDevice);
        bench_findOrInsert_kernel<<<numBlocks, blockSize>>>(view, d_keys, d_out, N);
        cudaMemcpy(h_out.data(), d_out, N * sizeof(uint), cudaMemcpyDeviceToHost);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&pipeline_ms, start, stop);

        cudaFree(d_keys);
        cudaFree(d_out);
        cudaFree(d_preload_keys);
        cudaFree(d_preload_out);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);

        return {kernel_ms, pipeline_ms};
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Run complete benchmark with configuration
        @param cfg Benchmark configuration
        @return Result structure with timing data */
    Result runBenchmarkInternal(const BenchmarkConfig& cfg)
    {
        Result R{};
        int    capacity  = static_cast<int>(cfg.capacity);
        int    repeats   = static_cast<int>(cfg.numTrials);
        int    N_preload = static_cast<int>(std::max(
            1u,
            static_cast<unsigned>(std::floor(cfg.loadFactor * static_cast<double>(capacity)))));
        int    N         = static_cast<int>(cfg.numLookups);

        if(N <= 0)
            N = 1;
        if(N_preload <= 0)
            N_preload = 1;
        R.N = N;

        // Generate keys for pre-populating the table
        auto preload_keys = makeKeys(N_preload, cfg.randomSeed);
        // Generate keys for lookup benchmarking (different seed to avoid exact matches)
        auto keys = makeKeys(N, cfg.randomSeed + 1);

        // Host: pre-populate table, then benchmark lookups (average over repeats)
        double host_total = 0.0;
        for(int r = 0; r < repeats; ++r)
        {
            host_total += runHostOnce(N, keys, N_preload, preload_keys, capacity);
        }
        R.host_ms = host_total / static_cast<double>(repeats);

        // Device: pre-populate table, then benchmark lookups (average over repeats)
        double kernel_total   = 0.0;
        double pipeline_total = 0.0;
        for(int r = 0; r < repeats; ++r)
        {
            auto [k_ms, p_ms] = runDeviceOnce(N, keys, N_preload, preload_keys, capacity);
            kernel_total += k_ms;
            pipeline_total += p_ms;
        }
        R.kernel_ms   = kernel_total / static_cast<double>(repeats);
        R.pipeline_ms = pipeline_total / static_cast<double>(repeats);

        return R;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Append result to CSV file
        @param result Result structure to write */
    void appendCSV(const Result& result)
    {
        // Write config fields followed by timing results for consistency
        m_csvWriter->writeRow(m_config.numTrials,
                              m_config.capacity,
                              m_config.loadFactor,
                              m_config.numLookups,
                              result.host_ms,
                              result.kernel_ms,
                              result.pipeline_ms);
    }
};

#endif