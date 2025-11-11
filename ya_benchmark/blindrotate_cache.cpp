//
// Created by xsong93 on 2025/11/11.
//

#include <nlohmann/json.hpp>

#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/cache_manager.h"
#include "yautil/cache_workload_generator.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

using json = nlohmann::json;

void benchStat(const std::vector<long>& iteration_times_us, const string& benchName, const string& saveFileName) {

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
    nlohmann::json results;
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

void benchLazy(const YatfheParameters& param, SimpleCacheManager& cache, const vector<int>& accessPattern) {
    // client side
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
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    // data gen
    Integer pt = 3;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);


    std::vector<long> iterationTimesUs;
    for (auto i = 1; i <= 50; i++) {
        clearFileCache();
        auto id = accessPattern[i];
        std::string file = DiskReader::generateLazyKeyFilename(id);
        auto start = steady_clock::now();

        ScaledTlwe sTlwe {param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        auto* bskServer = cache.getLazyKeySimple(id);
        Trlwe out{param};
        if (bskServer != nullptr) {
            blindRotateLazyPipeAltNtt(out, bskServer->bskFirst, bskServer->bskPrime,bskServer->s2Dft, sTlwe, v, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            printMsg(file, "Loading Lazy key from disk");
            blindRotateLazyPipeAltInitNtt(out, bsk.bskFirst, bsk.bskPrime,bsk.s2Dft, sTlwe,
                v, file, param);
            cache.putLazyKey(id, bsk);
        }

        Tlwe tmp{ksKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, out, param.driftPhase);
        switchKeyForTlwe(output, ksKey, tmp, param);

        auto elapsedUs = duration_cast<microseconds>(steady_clock::now() - start).count();

        iterationTimesUs.push_back(elapsedUs);
    }
    benchStat(iterationTimesUs, "Benchmark/LAZY", "lazy_benchmark_results.json");
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.lazyRequest << std::endl;
    std::cout << "Hit rate: " << stats.lazyHitRate() * 100 << "%" << std::endl;
}

void benchGinx(const YatfheParameters& param, SimpleCacheManager& cache, const vector<int>& accessPattern) {
    // client side
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
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    // data gen
    Integer pt = 3;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);

    // server side
    std::vector<long> iterationTimesUs;
    for (auto i = 1; i <= 50; i++) {
        clearFileCache();
        auto start = steady_clock::now();

        ScaledTlwe sTlwe {param.N * 2, param.n};
        Trlwe acc{param};
        rescaleTlweToNewMod(sTlwe, input);
        genNoiselessTrlweSample(acc, v, sTlwe);
        auto& bskServer = cache.getGinxKey(accessPattern[i]);
        blindRotateJP22Ntt(acc, bskServer.bskDft, sTlwe, param);
        Tlwe tmp{ksKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, acc, param.driftPhase);
        switchKeyForTlwe(output, ksKey, tmp, param);

        auto elapsedUs = duration_cast<microseconds>(steady_clock::now() - start).count();

        iterationTimesUs.push_back(elapsedUs);
    }
    benchStat(iterationTimesUs, "Benchmark/GINX", "ginx_benchmark_results.json");
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.ginxRequest << std::endl;
    std::cout << "Hit rate: " << stats.ginxHitRate() * 100 << "%" << std::endl;
}

int main(int argc, char **argv) {
    YatfheParameters param{};
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    // default
    int cacheCapacity = 5;
    double zipfParam = 0.8;
    int patternSize = 1000;

    if (argc > 1) {
        cacheCapacity = std::atoi(argv[1]);
    }
    if (argc > 2) {
        zipfParam = std::atof(argv[2]);
    }
    if (argc > 3) {
        patternSize = std::atoi(argv[3]);
    }

    printf("User param: Cache capacity=%d, Zipf s=%.1f, Max request count=%d\n",
           cacheCapacity, zipfParam, patternSize);

    // init cache
    SimpleCacheManager cache(param.n, cacheCapacity);
    CacheWorkloadGenerator workload(zipfParam);
    auto accessPattern = workload.generateAccessPattern(patternSize);
    printArray(accessPattern, "access pattern");

    benchGinx(param, cache, accessPattern);
    benchLazy(param, cache, accessPattern);

    return 0;
}
