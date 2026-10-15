#include <iostream>
#include <sys/stat.h>

#include "BenchmarkConfig.hh"
#include "ContactTableBenchmark.hh"

int main(int argc, char** argv)
{
    // Parse minimal args into BenchmarkConfig (capacity, loadFactor, numLookups, numTrials)
    BenchmarkConfig cfg;
    if(argc > 1)
        cfg.capacity = static_cast<uint>(std::atoi(argv[1]));
    if(argc > 2)
        cfg.loadFactor = std::atof(argv[2]);
    if(argc > 3)
        cfg.numLookups = static_cast<uint>(std::atoi(argv[3]));
    if(argc > 4)
        cfg.numTrials = static_cast<uint>(std::atoi(argv[4]));
    std::string csv
        = (argc > 5) ? std::string(argv[5]) : std::string("data/bench_contact_table.csv");

    std::cout << "Running contact table benchmark capacity=" << cfg.capacity
              << " loadFactor=" << cfg.loadFactor << " numLookups=" << cfg.numLookups
              << " numTrials=" << cfg.numTrials << "\n";

    // Determine if this is the first run (create new CSV) or append mode
    struct stat st;
    bool        appendMode = (stat(csv.c_str(), &st) == 0);

    // Create benchmark instance and run
    ContactTableBenchmark benchmark(cfg, csv, appendMode);
    benchmark.runBenchmark();

    std::cout << "Results appended to " << csv << "\n";

    return 0;
}
