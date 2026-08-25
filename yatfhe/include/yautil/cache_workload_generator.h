//
// Created by xsong93 on 11/10/25.
//

#ifndef YATFHE_CACHE_WORKLOAD_GENERATOR_H
#define YATFHE_CACHE_WORKLOAD_GENERATOR_H

#include <cstdint>
#include <vector>
#include <random>
#include <algorithm>
#include "yautil/zipf.h"

class CacheWorkloadGenerator {
private:
    ZipfDistribution zipf;
    int idBase;

public:
    explicit CacheWorkloadGenerator(double zipfParam = 1.2, int dataSize = 50,
                                    uint64_t seed = ZipfDistribution::kDefaultSeed, int idBase = 1)
        : zipf(zipfParam, dataSize, seed), idBase(idBase) {}

    uint64_t seed() const { return zipf.seed(); }

    int userCount() const { return zipf.size(); }

    // generate access pattern with Zipf distribution
    std::vector<int> generateAccessPattern(size_t length) {
        std::vector<int> pattern = zipf.generateBatch(length);
        for (int& id : pattern) {
            id += idBase;
        }
        return pattern;
    }

    // hot, cold data
    std::vector<int> generate_realistic_workload(size_t total_accesses, size_t dataSize, double hot_ratio = 0.8) {
        size_t hot_accesses = total_accesses * hot_ratio;
        size_t cold_accesses = total_accesses - hot_accesses;

        std::vector<int> workload;
        workload.reserve(total_accesses);

        // hot data (80% access)
        int hot_range = dataSize * 0.2;
        std::uniform_int_distribution<int> hot_dist(0, hot_range - 1);

        for (size_t i = 0; i < hot_accesses; ++i) {
            workload.push_back(hot_dist(rng));
        }

        // cold data（20% access）
        std::uniform_int_distribution<int> cold_dist(hot_range, dataSize);

        for (size_t i = 0; i < cold_accesses; ++i) {
            workload.push_back(cold_dist(rng));
        }

        std::shuffle(workload.begin(), workload.end(), rng);
        return workload;
    }
};

#endif //YATFHE_CACHE_WORKLOAD_GENERATOR_H
