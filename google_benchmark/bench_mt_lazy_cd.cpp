#include <thread>
#include <benchmark/benchmark.h>
#include <vector>
#include <utility>
#include <iostream>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"

class BlindRotateBenchmark : public benchmark::Fixture {
public:
    BlindRotateBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        param = YatfheParameters{};
        param.N = 512;
        initYatfhe(param);

        // key gen
        TlweKey tlweKey{param.n, param.lweStdDev};
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        TlweKeySwitchingKey ksKey{param};
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
        s2Dft = TrlevDft{param, param.l};
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        bskMPLazy = BootstrappingKeyMPOpt{param, param.lApprox, true};
        genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        sTlwe = ScaledTlwe{param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        acc = Trlwe{param.k, param.N};
        genNoiselessTrlweSample(acc, v, sTlwe);

        out = Trlwe{param.k, param.N};
    }

protected:
    YatfheParameters param;
    BootstrappingKeyMPOpt bskMPLazy;
    Trlwe acc;
    Trlwe out;
    TrlevDft s2Dft;
    ScaledTlwe sTlwe;
    TorusPolynomial v;
};

// Standard benchmark for comparison
BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_MULTITHREAD)(benchmark::State& state) {
    const int batchSize = state.range(0);
    const int tasksPerThread = state.range(1);
    auto localParam = param;
    localParam.batchSize = batchSize;
    localParam.tasksPerThread = tasksPerThread;

    for (auto _ : state) {
        blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft, bskMPLazy.initialized, sTlwe, v, s2Dft, localParam);
        benchmark::DoNotOptimize(out);
    }

    state.counters["batchSize"] = batchSize;
    state.counters["tasksPerThread"] = tasksPerThread;
}

// Coordinate descent optimized benchmark
BENCHMARK_DEFINE_F(BlindRotateBenchmark, OPTIMIZED_MULTITHREAD)(benchmark::State& state) {
    // Initial parameters
    int batchSize = 1;
    int tasksPerThread = 1;
    double bestThroughput = 0.0;

    // Coordinate descent parameters
    const int maxIterations = 5;
    const double improvementThreshold = 0.01; // 1% improvement

    // We'll store all benchmark runs for analysis
    std::vector<std::pair<std::pair<int, int>, double>> all_runs;

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        // Optimize tasksPerThread with fixed batchSize
        double bestTptThroughput = 0.0;
        int bestTpt = tasksPerThread;

        for (int tpt = 1; tpt <= 20; ++tpt) {
            auto localParam = param;
            localParam.batchSize = batchSize;
            localParam.tasksPerThread = tpt;

            // Measure performance
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 10; i++) {
                blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft,
                                  bskMPLazy.initialized, sTlwe, v, s2Dft, localParam);
                benchmark::DoNotOptimize(out);
            }
            auto end = std::chrono::high_resolution_clock::now();

            double duration = std::chrono::duration<double>(end - start).count();
            double throughput = 10.0 / duration; // operations per second

            all_runs.emplace_back(std::make_pair(batchSize, tpt), throughput);

            if (throughput > bestTptThroughput) {
                bestTptThroughput = throughput;
                bestTpt = tpt;
            }
        }

        tasksPerThread = bestTpt;

        // Optimize batchSize with fixed tasksPerThread
        double best_bs_throughput = 0.0;
        int best_bs = batchSize;

        for (int bs = 1; bs <= 20; ++bs) {
            auto localParam = param;
            localParam.batchSize = bs;
            localParam.tasksPerThread = tasksPerThread;

            // Measure performance
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 10; i++) {
                blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft,
                                  bskMPLazy.initialized, sTlwe, v, s2Dft, localParam);
                benchmark::DoNotOptimize(out);
            }
            auto end = std::chrono::high_resolution_clock::now();

            double duration = std::chrono::duration<double>(end - start).count();
            double throughput = 10.0 / duration; // operations per second

            all_runs.emplace_back(std::make_pair(bs, tasksPerThread), throughput);

            if (throughput > best_bs_throughput) {
                best_bs_throughput = throughput;
                best_bs = bs;
            }
        }

        // Check for convergence
        if (best_bs_throughput - bestThroughput < improvementThreshold * bestThroughput) {
            break;
        }

        batchSize = best_bs;
        bestThroughput = best_bs_throughput;
    }

    // Run final optimized configuration for the benchmark
    auto localParam = param;
    localParam.batchSize = batchSize;
    localParam.tasksPerThread = tasksPerThread;

    for (auto _ : state) {
        blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft,
                          bskMPLazy.initialized, sTlwe, v, s2Dft, localParam);
        benchmark::DoNotOptimize(out);
    }

    // Output optimization results
    std::cout << "\nOptimization Results:\n";
    std::cout << "Final batchSize: " << batchSize << "\n";
    std::cout << "Final tasksPerThread: " << tasksPerThread << "\n";
    std::cout << "Best throughput: " << bestThroughput << " ops/s\n";

    std::cout << "\nAll Tested Configurations:\n";
    for (const auto& run : all_runs) {
        std::cout << "batchSize=" << run.first.first
                  << ", tasksPerThread=" << run.first.second
                  << ": " << run.second << " ops/s\n";
    }

    state.counters["batchSize"] = batchSize;
    state.counters["tasksPerThread"] = tasksPerThread;
    state.counters["throughput"] = benchmark::Counter(bestThroughput,
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
}

// Register benchmarks
BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->ArgsProduct({
        benchmark::CreateDenseRange(1, 20, 1),  // batchSize
        benchmark::CreateDenseRange(1, 20, 1)   // tasksPerThread
    })
    ->ArgNames({"batchSize", "tasksPerThread"})
    ->MeasureProcessCPUTime()
    ->UseRealTime();

BENCHMARK_REGISTER_F(BlindRotateBenchmark, OPTIMIZED_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->MeasureProcessCPUTime()
    ->UseRealTime();

BENCHMARK_MAIN();