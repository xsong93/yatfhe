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

class CommandLineParser {
private:
    std::map<std::string, std::string> arguments_;
    std::vector<std::string> positionalArgs_;

public:
    CommandLineParser(int argc, char** argv) {
        parseArguments(argc, argv);
    }

    int getInt(const std::string& name, int defaultValue) const {
        auto it = arguments_.find(name);
        if (it != arguments_.end()) {
            return std::atoi(it->second.c_str());
        }
        return defaultValue;
    }

    double getDouble(const std::string& name, double defaultValue) const {
        auto it = arguments_.find(name);
        if (it != arguments_.end()) {
            return std::atof(it->second.c_str());
        }
        return defaultValue;
    }

    std::string getString(const std::string& name, const std::string& defaultValue) const {
        auto it = arguments_.find(name);
        if (it != arguments_.end()) {
            return it->second;
        }
        return defaultValue;
    }

    bool hasFlag(const std::string& flag) const {
        return arguments_.find(flag) != arguments_.end();
    }

    void printUsage(const std::string& programName) const {
        std::cout << "Usage: " << programName << " [Options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  --cap=N      Cache capacity (Default: 5)" << std::endl;
        std::cout << "  --s=N        Zipf s (Default: 0.8)" << std::endl;
        std::cout << "  --m=N        Size ratio of two keys (Default: 4)" << std::endl;
        std::cout << "  --pat=N      Request size (Default: 1000)" << std::endl;
        std::cout << "  --help       Show helps." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage example:" << std::endl;
        std::cout << "  " << programName << " --cap=10 --s=1.0 --pat=2000" << std::endl;
        std::cout << "  " << programName << " --s=0.9" << std::endl;
    }

private:
    void parseArguments(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg == "--help") {
                arguments_["help"] = "true";
            }
            else if (arg.substr(0, 2) == "--") {
                size_t pos = arg.find('=');
                if (pos != std::string::npos) {
                    std::string name = arg.substr(2, pos - 2);
                    std::string value = arg.substr(pos + 1);
                    arguments_[name] = value;
                } else {
                    arguments_[arg.substr(2)] = "true";
                }
            }
            else {
                positionalArgs_.push_back(arg);
            }
        }
    }
};

void benchStat(const std::vector<long>& iteration_times_us, const long request, const double hitRate, const string& benchName,
               const int cacheCap, const string& saveFileName) {

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
    results["pressure_ratio"] = 50.0 / cacheCap;
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
        {"total_time_us", sum},
        {"total_requests", request},
        {"hit_rate", hitRate}
    };

    // Write to file
    std::ofstream outfile(saveFileName);
    outfile << results.dump(4) << std::endl; // Pretty print with 4-space indent
    outfile.close();

    std::cout << "Benchmark results saved to: " << saveFileName << std::endl;
    std::cout << "Average time: " << average_time << " μs" << std::endl;
}

void benchLazy(const YatfheParameters& param, SimpleCacheManager& cache, const vector<int>& accessPattern,
               const int cacheCap, bool isSave) {
    cout << "bench lazy" << endl;
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


    // server
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe out{param};

    // warm up
    cout << "warm up" << endl;
    for (auto i = 0; i < accessPattern.size()/2; i++) {
        cache.getLazyKey(accessPattern[i]);
        steady_clock::now();
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<long> iterationTimesUs;
    for (auto i = accessPattern.size()/2; i < 100 + accessPattern.size()/2; i++) {
        clearFileCache();
        auto id = accessPattern[i];
        std::string file = DiskReader::generateLazyKeyFilename(id);
        auto start = steady_clock::now();
        auto* bskServer = cache.getLazyKeySimple(id);
        if (bskServer != nullptr) {
            blindRotateLazyPipeAltNtt(out, bskServer->bskFirst, bskServer->bskPrime,bskServer->s2Dft, sTlwe, v, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            blindRotateLazyPipeAltInitNtt(out, bsk.bskFirst, bsk.bskPrime,bsk.s2Dft, sTlwe, v, file, param);
            cache.putLazyKey(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        auto elapsedUs = duration_cast<microseconds>(end - start).count();
        iterationTimesUs.push_back(elapsedUs);
    }
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.lazyRequest << std::endl;
    std::cout << "Hit rate: " << stats.lazyHitRate() * 100 << "%" << std::endl;
    string file = "ours_benchmark_results_" + to_string(cacheCap) + ".json";
    if (isSave) {
        benchStat(iterationTimesUs, iterationTimesUs.size(), stats.lazyHitRate(),
                  "Benchmark/OURS",  cacheCap, file);
    }
}

void benchGinx(const YatfheParameters& param, SimpleCacheManager& cache, const vector<int>& accessPattern,
               const int cacheCap, bool isSave) {
    cout << "bench ginx" << endl;
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
    ScaledTlwe sTlwe {param.N * 2, param.n};
    Trlwe acc{param};
    rescaleTlweToNewMod(sTlwe, input);
    genNoiselessTrlweSample(acc, v, sTlwe);

    // warm up
    cout << "warm up" << endl;
    for (auto i = 0; i < accessPattern.size()/2; i++) {
        cache.getGinxKey(accessPattern[i]);
        steady_clock::now();
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<long> iterationTimesUs;
    for (auto i = accessPattern.size()/2; i < 100 + accessPattern.size()/2; i++) {
        clearFileCache();
        auto id = accessPattern[i];
        std::string file = DiskReader::generateGinxKeyFilename(id);
        auto start = steady_clock::now();
        // auto& bskServer = cache.getGinxKey(accessPattern[i]);
        // blindRotateJP22Ntt(acc, bskServer.bskDft, sTlwe, param);
        auto* bskServer = cache.getGinxKeySimple(id);
        if (bskServer != nullptr) {
            blindRotateJP22Ntt(acc, bskServer->bskDft, sTlwe, param);
        } else {
            BootstrappingKeyMP bsk;
            deserializeBskMP(bsk, file, param.n);
            blindRotateJP22Ntt(acc, bsk.bskDft, sTlwe, param);
            cache.putMpKey(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        auto elapsedUs = duration_cast<microseconds>(end - start).count();
        iterationTimesUs.push_back(elapsedUs);
    }
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.ginxRequest << std::endl;
    std::cout << "Hit rate: " << stats.ginxHitRate() * 100 << "%" << std::endl;
    string file = "tfhe_benchmark_results_" + to_string(cacheCap) + ".json";
    if (isSave) {
        benchStat(iterationTimesUs, iterationTimesUs.size(), stats.ginxHitRate(),
                  "Benchmark/TFHE", cacheCap, file);
    }
}

void benchWWL24(const YatfheParameters& param, SimpleCacheManager& cache, const vector<int>& accessPattern,
                const int cacheCap, bool isSave) {
    cout << "bench WWL24" << endl;
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
    ScaledTlwe sTlwe {param.N * 2, param.n};
    Trlwe acc{param};
    rescaleTlweToNewMod(sTlwe, input);
    genNoiselessTrlweSample(acc, v, sTlwe);

    // warm up
    cout << "warm up" << endl;
    for (auto i = 0; i < accessPattern.size()/2; i++) {
        cache.getWWL24Key(accessPattern[i]);
        steady_clock::now();
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<long> iterationTimesUs;
    for (auto i = accessPattern.size()/2; i < 100 + accessPattern.size()/2; i++) {
        clearFileCache();
        auto id = accessPattern[i];
        std::string file = DiskReader::generateWWL24KeyFilename(id);
        auto start = steady_clock::now();
        auto* bskServer = cache.getWWL24KeySimple(id);
        if (bskServer != nullptr) {
            blindRotateWWL24Ntt(acc, bskServer->bskDft, sTlwe, bskServer->s2Dft, param);
        } else {
            BootstrappingKeyWWL24 bsk;
            deserializeBskWWL24(bsk, file, param.n);
            blindRotateWWL24Ntt(acc, bsk.bskDft, sTlwe, bsk.s2Dft, param);
            cache.putWWL24Key(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        auto elapsedUs = duration_cast<microseconds>(end - start).count();
        iterationTimesUs.push_back(elapsedUs);
    }
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.wwl24Request << std::endl;
    std::cout << "Hit rate: " << stats.wwl24HitRate() * 100 << "%" << std::endl;
    string file = "wwl+24_benchmark_results_" + to_string(cacheCap) + ".json";
    if (isSave) {
        benchStat(iterationTimesUs, iterationTimesUs.size(), stats.wwl24HitRate(),
                  "Benchmark/WWL+24", cacheCap, file);
    }
}

int main(int argc, char **argv) {
    int sizeRatio = round(67156489.0/16919095.0); // ginx key size / lazy key size
    int sizeRatio2 = round(67156489.0/33709652.0); // ginx key size / wwl+24 key size

    CommandLineParser parser(argc, argv);

    if (parser.hasFlag("help")) {
        parser.printUsage(argv[0]);
        return 0;
    }

    int cacheCapacity = parser.getInt("cap", 5);
    double zipfParam = parser.getDouble("s", 0.83);
    int multiplier = parser.getInt("m", sizeRatio);
    int multiplier2 = parser.getInt("m2", sizeRatio2);
    int patternSize = parser.getInt("pat", 2000);

    if (cacheCapacity <= 0) {
        std::cerr << "Error: Cache capacity should be larger than 0，using default 5" << std::endl;
        cacheCapacity = 5;
    }
    if (zipfParam <= 0) {
        std::cerr << "Error: Zipf parameter should be larger than 0，using default 0.83" << std::endl;
        zipfParam = 0.83;
    }
    if (multiplier <= 0) {
        std::cerr << "Error: m should be no less than 1，using default 4" << std::endl;
        multiplier = sizeRatio;
    }
    if (patternSize <= 0) {
        std::cerr << "Error: Request pattern size should be larger than 0，using default 2000" << std::endl;
        patternSize = 2000;
    }

    // printf("User param: Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n",
    //        cacheCapacity, zipfParam, patternSize);

    YatfheParameters param{};
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);


    // warm up cycle
    {
        cacheCapacity = 100;
        multiplier = 1;
        multiplier2 = 1;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    // benchmarking
    {
        cacheCapacity = 50;
        multiplier = 1;
        multiplier2 = 1;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    {
        cacheCapacity = 25;
        multiplier = sizeRatio;
        multiplier2 = sizeRatio2;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    {
        cacheCapacity = 10;
        multiplier = sizeRatio;
        multiplier2 = sizeRatio2;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    {
        cacheCapacity = 5;
        multiplier = sizeRatio;
        multiplier2 = sizeRatio2;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    {
        cacheCapacity = 1;
        multiplier = 1;
        multiplier2 = 1;
        printf("Cache capacity=%d, Zipf s=%.3f, Max request count=%d\n", cacheCapacity, zipfParam, patternSize);
        // init cache
        SimpleCacheManager cache(param.n, cacheCapacity, cacheCapacity * multiplier, cacheCapacity * multiplier2);
        CacheWorkloadGenerator workload(zipfParam);
        auto accessPattern = workload.generateAccessPattern(patternSize);
        printArray(accessPattern, "access pattern");

        benchWWL24(param, cache, accessPattern, cacheCapacity, true);
        benchLazy(param, cache, accessPattern, cacheCapacity, true);
        benchGinx(param, cache, accessPattern, cacheCapacity, true);
    }

    return 0;
}
