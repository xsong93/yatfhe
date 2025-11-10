#include <chrono>
#include <benchmark/benchmark.h>
#include <nlohmann/json.hpp>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/cache_manager.h"
#include "yautil/cache_workload_generator.h"
#include "yautil/lru_cache.h"
#include "yautil/tool.h"

using json = nlohmann::json;

class ZipfBenchmark : public benchmark::Fixture {
public:
    ZipfBenchmark()
        : param{},
        v{param.N},
        tlweKey{param.n, param.lweStdDev},
        trgswKey{param},
        dummyKsKey{param},
        input{param.n},
        acc{param},
        sTlwe{param.N * 2, param.n},
        cache{param.n,10},
        workload{1.0},
        accessPattern{workload.generateAccessPattern(1000)} {}

    void SetUp(const benchmark::State& state) override {
        initYatfhe(param);

        // key gen
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        generateTestPolynomial(v, param.torusBase, 2 * param.N);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(dummyKsKey, trlweKey, tlweKsKey, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        symEncTlwe(input, mu, tlweKey);
    }

protected:
    YatfheParameters param;
    TorusPolynomial v;
    TlweKey tlweKey;
    TrgswKey trgswKey;
    TlweKeySwitchingKey dummyKsKey;
    Tlwe input;
    Trlwe acc;
    ScaledTlwe sTlwe;
    SimpleCacheManager cache;
    CacheWorkloadGenerator workload;
    vector<int> accessPattern;
};

void benchStat(const std::vector<double>& iteration_times_us, const string& benchName, const string& saveFileName) {

    // Calculate statistics
    double sum = 0.0;
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::min();

    for (double time : iteration_times_us) {
        sum += time;
        if (time < min_time) min_time = time;
        if (time > max_time) max_time = time;
    }

    double average_time = sum / iteration_times_us.size();

    // Create JSON structure
    json results;
    results["benchmark_name"] = benchName;
    results["total_iterations"] = iteration_times_us.size();
    results["time_unit"] = "microseconds";

    // Individual iteration times
    results["iterations"] = json::array();
    for (size_t i = 0; i < iteration_times_us.size(); ++i) {
        results["iterations"].push_back({
            {"iteration", i + 1},
            {"time_us", iteration_times_us[i]}
        });
    }

    // Statistics
    results["statistics"] = {
        {"average_time_us", average_time},
        {"min_time_us", min_time},
        {"max_time_us", max_time},
        {"total_time_us", sum}
    };

    // Write to file
    std::ofstream outfile(saveFileName);
    outfile << results.dump(4) << std::endl; // Pretty print with 4-space indent
    outfile.close();

    std::cout << "Benchmark results saved to: " << saveFileName << std::endl;
    std::cout << "Average time: " << average_time << " μs" << std::endl;
}

BENCHMARK_DEFINE_F(ZipfBenchmark, GINX)(benchmark::State& state) {
    std::vector<double> iterationTimesUs;

    for (auto _ : state) {
        auto start = std::chrono::high_resolution_clock::now();

        rescaleTlweToNewMod(sTlwe, input);
        genNoiselessTrlweSample(acc, v, sTlwe);
        auto& bskServer = cache.getGinxKey(accessPattern[state.iterations()]);
        blindRotateJP22Ntt(acc, bskServer.bskDft, sTlwe, param);
        Tlwe tmp{dummyKsKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, acc, param.driftPhase);
        switchKeyForTlwe(output, dummyKsKey, tmp, param);

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        iterationTimesUs.push_back(static_cast<double>(elapsedUs));

        benchmark::DoNotOptimize(output);
    }
    benchStat(iterationTimesUs, "Benchmark/GINX", "ginx_benchmark_results.json");
}

BENCHMARK_DEFINE_F(ZipfBenchmark, LAZY)(benchmark::State& state) {
    std::vector<double> iterationTimesUs;
    for (auto _ : state) {
        auto start = std::chrono::high_resolution_clock::now();

        rescaleTlweToNewMod(sTlwe, input);
        auto* bskServer = cache.getLazyKey(accessPattern[state.iterations()]);
        Trlwe out{param};
        if (bskServer != nullptr) {
            blindRotateLazyPipeAltNtt(out, bskServer->bskFirst, bskServer->bskPrime,bskServer->s2Dft, sTlwe, v, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            std::string file = DiskReader::generateLazyKeyFilename(state.iterations());
            blindRotateLazyPipeAltInitNtt(out, bsk.bskFirst, bsk.bskPrime,bsk.s2Dft, sTlwe,
                v, file, param);
            cache.putLazyKey(accessPattern[state.iterations()], bsk);
        }

        Tlwe tmp{dummyKsKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, out, param.driftPhase);
        switchKeyForTlwe(output, dummyKsKey, tmp, param);

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        iterationTimesUs.push_back(static_cast<double>(elapsedUs));

        benchmark::DoNotOptimize(output);
    }
    benchStat(iterationTimesUs, "Benchmark/LAZY", "lazy_benchmark_results.json");
}

BENCHMARK_REGISTER_F(ZipfBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseManualTime();
BENCHMARK_REGISTER_F(ZipfBenchmark, LAZY)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseManualTime();

BENCHMARK_MAIN();