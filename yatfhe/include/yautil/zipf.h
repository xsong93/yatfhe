//
// Created for anonymous review.
//

#ifndef BASE_ZIPF_H
#define BASE_ZIPF_H

#include <cstdint>
#include <vector>
#include <random>

class ZipfDistribution {
private:
    int N;
    std::vector<double> probabilities;
    std::vector<double> cumulative;
    std::uniform_real_distribution<double> dist;
    std::random_device rd_;
    std::mt19937 gen_;

public:
    ZipfDistribution(double s = 1.0, int N = 50) : gen_(rd_()), dist(0.0, 1.0), N(N) {
        probabilities.resize(N);
        cumulative.resize(N);
        calculateDistribution(s);
    }

    // generate random number
    int generate() {
        double u = dist(gen_);

        for (int i = 1; i <= N; ++i) {
            if (u <= cumulative[i]) {
                return i;
            }
        }
        return N - 1;
    }

    // get probability for a number
    double probability(int number) const {
        if (number < 0 || number >= N) return 0.0;
        return probabilities[number];
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
            int number = rank - 1;
            probabilities[number] = 1.0 / (std::pow(rank, s) * harmonic);
        }

        cumulative[0] = probabilities[0];
        for (int i = 1; i < N; ++i) {
            cumulative[i] = cumulative[i - 1] + probabilities[i];
        }

        cumulative[N - 1] = 1.0;
    }
};

#endif //BASE_ZIPF_H
