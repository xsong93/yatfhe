//
// Created by xsong93 on 11/10/25.
//

#ifndef YATFHE_ZIPF_H
#define YATFHE_ZIPF_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

class ZipfDistribution {
private:
    int N;
    uint64_t seed_;
    std::vector<double> probabilities;
    std::vector<double> cumulative;
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 gen_;

public:
    static constexpr uint64_t kDefaultSeed = 20260816ULL;

    explicit ZipfDistribution(double s = 1.0, int N = 50, uint64_t seed = kDefaultSeed)
        : N(N), seed_(seed), dist(0.0, 1.0), gen_(seed) {
        probabilities.resize(N);
        cumulative.resize(N);
        calculateDistribution(s);
    }

    static uint64_t randomSeed() {
        std::random_device rd;
        return (static_cast<uint64_t>(rd()) << 32) ^ rd();
    }

    uint64_t seed() const { return seed_; }

    int size() const { return N; }

    void reseed(const uint64_t seed) {
        seed_ = seed;
        gen_.seed(seed);
    }

    int generate() {
        const double u = dist(gen_);
        const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), u);
        auto rank = static_cast<int>(it - cumulative.begin());
        if (rank >= N) rank = N - 1;   // guards u == 1.0 and fp rounding at the tail
        return rank;
    }

    double probability(const int rank) const {
        if (rank < 0 || rank >= N) return 0.0;
        return probabilities[rank];
    }

    // generate random distributed access pattern
    std::vector<int> generateBatch(size_t count) {
        std::vector<int> result;
        result.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            result.push_back(generate());
        }
        return result;
    }

private:
    void calculateDistribution(double s) {
        double harmonic = 0.0;
        for (int k = 1; k <= N; ++k) {
            harmonic += 1.0 / std::pow(k, s);
        }

        for (int rank = 1; rank <= N; ++rank) {
            probabilities[rank - 1] = 1.0 / (std::pow(rank, s) * harmonic);
        }

        cumulative[0] = probabilities[0];
        for (int i = 1; i < N; ++i) {
            cumulative[i] = cumulative[i - 1] + probabilities[i];
        }

        cumulative[N - 1] = 1.0;
    }
};

#endif //YATFHE_ZIPF_H
