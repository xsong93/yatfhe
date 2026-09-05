
#include "include/bench_out.h"
#include <nlohmann/json.hpp>

#include <sched.h>
#include <unistd.h>
#include <algorithm>
#include <numeric>

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
        std::cout << "  --caps=LIST  Comma separated cache capacities (Default: 100,50,20,10,5,1)" << std::endl;
        std::cout << "  --users=N    Number of tenants (Default: 50)" << std::endl;
        std::cout << "  --s=N        Zipf skew s (Default: 0.83)" << std::endl;
        std::cout << "  --reqs=N     Measured requests per method (Default: 10000)" << std::endl;
        std::cout << "  --warm=N     Warm-up requests per method (Default: 500)" << std::endl;
        std::cout << "  --seed=N     Workload RNG seed (Default: 20260816)" << std::endl;
        std::cout << "  --iso=MODE   Cache budget model: bytes | slots (Default: bytes)" << std::endl;
        std::cout << "               bytes: equal memory" << std::endl;
        std::cout << "               slots: equal slot count" << std::endl;
        std::cout << "  --m=N        GINX/OURS key size ratio (Default: 4)" << std::endl;
        std::cout << "  --m2=N       GINX/WWL+24 key size ratio (Default: 2)" << std::endl;
        std::cout << "  --tag=STR    Suffix appended to the output json" << std::endl;
        std::cout << "  --out=DIR    Directory for result json" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage example:" << std::endl;
        std::cout << "  taskset -c 0-7 " << programName << " --reqs=10000 --seed=1" << std::endl;
        std::cout << "  " << programName << " --caps=10,5 --reqs=2000 --iso=slots --tag=iso" << std::endl;
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

struct RunConfig {
    int users{};
    int ginxCapacity{};
    int pipeCapacity{};
    int wwl24Capacity{};
    double zipfS{};
    uint64_t seed{};
    long warmupRequests{};
    long measuredRequests{};
    std::string isoMode;
    std::string tag;
};

struct RequestRecord {
    long timeUs;
    int keyId;
    bool hit;
};

static std::string cpuAffinity() {
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) != 0) {
        return "none";
    }
    std::string out;
    int tmp = -1;
    const int maxCpu = CPU_SETSIZE;
    for (int cpu = 0; cpu <= maxCpu; ++cpu) {
        const bool in = cpu < maxCpu && CPU_ISSET(cpu, &mask);
        if (in && tmp < 0) {
            tmp = cpu;
        } else if (!in && tmp >= 0) {
            if (!out.empty())
                out += ',';
            out += std::to_string(tmp);
            if (cpu - 1 != tmp)
                out += "-" + std::to_string(cpu - 1);
            tmp = -1;
        }
    }
    return out.empty() ? "none" : out;
}

static double percentile(const std::vector<long>& sorted, const double p) {
    if (sorted.empty())
        return 0.0;
    const double rank = p / 100.0 * static_cast<double>(sorted.size());
    auto idx = static_cast<long>(std::ceil(rank)) - 1;
    if (idx < 0)
        idx = 0;
    if (idx >= static_cast<long>(sorted.size()))
        idx = static_cast<long>(sorted.size()) - 1;
    return static_cast<double>(sorted[idx]);
}

static json summarize(std::vector<long> samples) {
    json out;
    out["count"] = samples.size();
    if (samples.empty()) {
        return out;
    }
    std::sort(samples.begin(), samples.end());
    const double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
    const double mean = sum / static_cast<double>(samples.size());
    double sq = 0.0;
    for (const long v : samples) {
        const double d = static_cast<double>(v) - mean;
        sq += d * d;
    }
    const double sd = std::sqrt(sq / static_cast<double>(samples.size()));

    out["mean_us"] = mean;
    out["sd_us"] = sd;
    out["min_us"] = static_cast<double>(samples.front());
    out["max_us"] = static_cast<double>(samples.back());
    out["p50_us"] = percentile(samples, 50.0);
    out["p90_us"] = percentile(samples, 90.0);
    out["p99_us"] = percentile(samples, 99.0);
    out["p99_9_us"] = percentile(samples, 99.9);
    return out;
}

void benchStat(const std::vector<RequestRecord>& records, const double hitRate,
               const string& benchName, const RunConfig& cfg, const string& saveFileName) {
    std::vector<long> all, hits, misses;
    all.reserve(records.size());
    for (const auto& [timeUs, keyId, hit] : records) {
        all.push_back(timeUs);
        (hit ? hits : misses).push_back(timeUs);
    }

    json results;
    results["benchmark_name"] = benchName;
    results["schema_version"] = 2;
    results["time_unit"] = "microseconds";

    results["config"] = {
        {"users", cfg.users},
        {"ginx_capacity_keys", cfg.ginxCapacity},
        {"lazy_capacity_keys", cfg.pipeCapacity},
        {"wwl24_capacity_keys", cfg.wwl24Capacity},
        {"pressure_ratio", static_cast<double>(cfg.users) / cfg.ginxCapacity},
        {"zipf_s", cfg.zipfS},
        {"seed", cfg.seed},
        {"warmup_requests", cfg.warmupRequests},
        {"measured_requests", cfg.measuredRequests},
        {"iso_mode", cfg.isoMode},
        {"tag", cfg.tag},
    };

    results["environment"] = {
        {"affinity_cpus", cpuAffinity()},
        {"cpus_available", static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN))},
    };

    results["iterations"] = json::array();
    for (size_t i = 0; i < records.size(); ++i) {
        results["iterations"].push_back({
            {"iteration", i + 1},
            {"time_us", records[i].timeUs},
            {"key_id", records[i].keyId},
            {"hit", records[i].hit},
        });
    }

    results["statistics"] = summarize(all);
    results["statistics"]["hit_rate"] = hitRate;
    results["statistics"]["total_requests"] = records.size();
    results["hit_statistics"] = summarize(hits);
    results["miss_statistics"] = summarize(misses);

    const std::string savePath = yabench::benchOutPath(saveFileName);
    std::ofstream outfile(savePath);
    outfile << results.dump(2) << std::endl;
    outfile.close();

    std::cout << "Saved " << savePath
              << "  n=" << all.size()
              << "  mean=" << results["statistics"]["mean_us"].get<double>() / 1000.0 << " ms"
              << "  sd=" << results["statistics"]["sd_us"].get<double>() / 1000.0 << " ms"
              << "  p99.9=" << results["statistics"]["p99_9_us"].get<double>() / 1000.0 << " ms"
              << "  hr=" << hitRate << std::endl;
}

static string outputName(const string& stem, const RunConfig& cfg) {
    string name = stem + "_benchmark_results_" + to_string(cfg.ginxCapacity);
    if (!cfg.tag.empty()) name += "_" + cfg.tag;
    return name + ".json";
}

void benchPipe(const Tlwe& input, const YatfheParameters& param, SimpleCacheManager& cache,
               const vector<int>& accessPattern, const RunConfig& cfg, bool isSave) {
    cout << "bench ours" << endl;

    // server
    TorusPolynomial v {param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);
    vector<NttPolynomial> gdVntt;
    prepareGdV(gdVntt, v, param);
    ScaledTlwe sTlwe {2 * param.N, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe out{param};

    // warm up: pre-touch every tenant once, then the Zipf warm-up sequence.
    cout << "warm up" << endl;
    for (int id = 1; id <= cfg.users; id++) {
        cache.getLazyKey(id);
    }
    for (long i = 0; i < cfg.warmupRequests; i++) {
        cache.getLazyKey(accessPattern[i]);
    }
    // settle the compute path before the measured window
    if (auto* key = cache.getLazyKeySimple(1)) {
        for (int k = 0; k < 20; k++) {
            blindRotateLazyPipeAltNtt(out, *key, sTlwe, v, gdVntt, param);
        }
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<RequestRecord> records;
    records.reserve(cfg.measuredRequests);
    for (long i = cfg.warmupRequests; i < cfg.warmupRequests + cfg.measuredRequests; i++) {
        auto id = accessPattern[i];
        std::string file = DiskReader::generateLazyKeyFilename(id);
        // Evict before the timer: on a cache miss this file is read.
        clearFileCache(file);
        auto start = steady_clock::now();
        auto* bskServer = cache.getLazyKeySimple(id);
        const bool hit = bskServer != nullptr;
        if (hit) {
            blindRotateLazyPipeAltNtt(out, *bskServer, sTlwe, v, gdVntt, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            blindRotateLazyPipeAltInitNtt(out, bsk, sTlwe, v, file, gdVntt, param);
            cache.putLazyKey(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        records.push_back({duration_cast<microseconds>(end - start).count(), id, hit});
    }
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.lazyRequest << std::endl;
    std::cout << "Hit rate: " << stats.lazyHitRate() * 100 << "%" << std::endl;
    if (isSave) {
        benchStat(records, stats.lazyHitRate(), "Benchmark/OURS", cfg, outputName("ours", cfg));
    }
}

void benchGinx(const Tlwe& input, const YatfheParameters& param, SimpleCacheManager& cache,
               const vector<int>& accessPattern, const RunConfig& cfg, bool isSave) {
    cout << "bench ginx" << endl;

    // server side
    ScaledTlwe sTlwe {2 * param.N, param.n};
    Trlwe acc{param};
    rescaleTlweToNewMod(sTlwe, input);
    TorusPolynomial v {param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    // warm up: pre-touch every tenant once, then the Zipf warm-up sequence.
    cout << "warm up" << endl;
    for (int id = 1; id <= cfg.users; id++) {
        cache.getGinxKey(id);
    }
    for (long i = 0; i < cfg.warmupRequests; i++) {
        cache.getGinxKey(accessPattern[i]);
    }
    // settle the compute path before the measured window
    if (auto* key = cache.getGinxKeySimple(1)) {
        for (int k = 0; k < 20; k++) {
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateJP22Ntt(acc, *key, sTlwe, param);
        }
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<RequestRecord> records;
    records.reserve(cfg.measuredRequests);
    for (long i = cfg.warmupRequests; i < cfg.warmupRequests + cfg.measuredRequests; i++) {
        auto id = accessPattern[i];
        std::string file = DiskReader::generateGinxKeyFilename(id);
        // Evict before the timer: on a cache miss this file is read.
        clearFileCache(file);
        genNoiselessTrlweSample(acc, v, sTlwe);
        auto start = steady_clock::now();
        auto* bskServer = cache.getGinxKeySimple(id);
        const bool hit = bskServer != nullptr;
        if (hit) {
            blindRotateJP22Ntt(acc, *bskServer, sTlwe, param);
        } else {
            BootstrappingKeyMP bsk;
            deserializeBskMP(bsk, file, param.n);
            blindRotateJP22Ntt(acc, bsk, sTlwe, param);
            cache.putMpKey(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        records.push_back({duration_cast<microseconds>(end - start).count(), id, hit});
    }
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.ginxRequest << std::endl;
    std::cout << "Hit rate: " << stats.ginxHitRate() * 100 << "%" << std::endl;
    if (isSave) {
        benchStat(records, stats.ginxHitRate(), "Benchmark/TFHE", cfg, outputName("tfhe", cfg));
    }
}

void benchWWL24(const Tlwe& input, const YatfheParameters& param, SimpleCacheManager& cache,
                const vector<int>& accessPattern, const RunConfig& cfg, bool isSave) {
    cout << "bench WWL24" << endl;

    // server side
    ScaledTlwe sTlwe {2 * param.N, param.n};
    Trlwe acc{param};
    rescaleTlweToNewMod(sTlwe, input);
    TorusPolynomial v {param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    // warm up: pre-touch every tenant once, then the Zipf warm-up sequence.
    cout << "warm up" << endl;
    for (int id = 1; id <= cfg.users; id++) {
        cache.getWWL24Key(id);
    }
    for (long i = 0; i < cfg.warmupRequests; i++) {
        cache.getWWL24Key(accessPattern[i]);
    }
    // settle the compute path before the measured window
    if (auto* key = cache.getWWL24KeySimple(1)) {
        for (int k = 0; k < 20; k++) {
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateWWL24Ntt(acc, *key, sTlwe, param);
        }
    }
    cache.resetStats();

    cout << "normal run" << endl;
    // normal run
    std::vector<RequestRecord> records;
    records.reserve(cfg.measuredRequests);
    for (long i = cfg.warmupRequests; i < cfg.warmupRequests + cfg.measuredRequests; i++) {
        auto id = accessPattern[i];
        std::string file = DiskReader::generateWWL24KeyFilename(id);
        clearFileCache(file);
        genNoiselessTrlweSample(acc, v, sTlwe);
        auto start = steady_clock::now();
        auto* bskServer = cache.getWWL24KeySimple(id);
        const bool hit = bskServer != nullptr;
        if (hit) {
            blindRotateWWL24Ntt(acc, *bskServer, sTlwe, param);
        } else {
            BootstrappingKeyWWL24 bsk;
            deserializeBskWWL24(bsk, file, param.n);
            blindRotateWWL24Ntt(acc, bsk, sTlwe, param);
            cache.putWWL24Key(id, std::move(bsk));
        }
        auto end = steady_clock::now();
        records.push_back({duration_cast<microseconds>(end - start).count(), id, hit});
    }
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.wwl24Request << std::endl;
    std::cout << "Hit rate: " << stats.wwl24HitRate() * 100 << "%" << std::endl;
    if (isSave) {
        benchStat(records, stats.wwl24HitRate(), "Benchmark/WWL+24", cfg, outputName("wwl+24", cfg));
    }
}

static std::vector<int> parseCapacity(const std::string& text) {
    std::vector<int> caps;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        const int value = std::atoi(item.c_str());
        if (value > 0) caps.push_back(value);
    }
    return caps;
}

int main(int argc, char **argv) {
    constexpr int sizeRatio = 4;  // ginx key size / ours key size
    constexpr int sizeRatio2 = 2; // ginx key size / wwl+24 key size

    CommandLineParser parser(argc, argv);

    if (parser.hasFlag("help")) {
        parser.printUsage(argv[0]);
        return 0;
    }

    const std::string caps = parser.getString("caps", "100,50,20,10,5,1");
    const int users = parser.getInt("users", 50);
    double zipfParam = parser.getDouble("s", 0.83);
    const int multiplier = parser.getInt("m", sizeRatio);
    const int multiplier2 = parser.getInt("m2", sizeRatio2);
    long requests = parser.getInt("reqs", 10000);
    long warmup = parser.getInt("warm", 500);
    const auto seed = static_cast<uint64_t>(parser.getInt("seed", ZipfDistribution::kDefaultSeed));
    const std::string iso = parser.getString("iso", "bytes");
    const std::string tag = parser.getString("tag", "");
    yabench::setBenchOutDir(parser.getString("out", yabench::benchOutDirRef()));

    std::vector<int> capacities = parseCapacity(caps);
    if (capacities.empty()) {
        std::cerr << "Error: --caps is empty, using the default" << std::endl;
        capacities = {100, 50, 20, 10, 5, 1};
    }
    if (zipfParam <= 0) {
        std::cerr << "Error: Zipf parameter should be larger than 0, using default 0.83" << std::endl;
        zipfParam = 0.83;
    }
    if (requests <= 0) {
        std::cerr << "Error: --reqs should be larger than 0, using default 10000" << std::endl;
        requests = 10000;
    }
    if (warmup < 0) {
        std::cerr << "Error: --warm should be larger than 0, using default 500" << std::endl;
        warmup = 500;
    }
    if (iso != "bytes" && iso != "slots") {
        std::cerr << "Error: --iso must be bytes or slots" << std::endl;
        return 1;
    }

    std::cout << "users=" << users
              << "zipf s=" << zipfParam
              << "seed=" << seed
              << "warm=" << warmup
              << "reqs=" << requests
              << "iso=" << iso << std::endl;

    YatfheParameters param{};
    initYatfhe(param);

    // client side
    // key gen
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    // data gen
    const int p = param.torusBase;
    const int encMod = 2 * p;
    Integer pt = -1;
    const int slot = static_cast<int>(((pt % p) + p) % p);
    Torus mu = modSwitchToTorusGeneral(slot, encMod, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);

    // server side
    const size_t patternSize = static_cast<size_t>(warmup + requests);
    for (const int cacheCapacity : capacities) {
        const int pipeMultiplier = iso == "bytes" ? multiplier : 1;
        const int wwl24Multiplier = iso == "bytes" ? multiplier2 : 1;

        RunConfig cfg;
        cfg.users = users;
        cfg.ginxCapacity = cacheCapacity;
        cfg.pipeCapacity = cacheCapacity * pipeMultiplier;
        cfg.wwl24Capacity = cacheCapacity * wwl24Multiplier;
        cfg.zipfS = zipfParam;
        cfg.seed = seed;
        cfg.warmupRequests = warmup;
        cfg.measuredRequests = requests;
        cfg.isoMode = iso;
        cfg.tag = tag;

        printf("\n=== capacity=%d (pressure %.2f)  ours=%d wwl24=%d slots  seed=%llu ===\n",
               cacheCapacity, static_cast<double>(users) / cacheCapacity,
               cfg.pipeCapacity, cfg.wwl24Capacity, static_cast<unsigned long long>(seed));

        SimpleCacheManager cache(param.n, cfg.ginxCapacity, cfg.pipeCapacity, cfg.wwl24Capacity);
        CacheWorkloadGenerator workload(zipfParam, users, seed);
        auto accessPattern = workload.generateAccessPattern(patternSize);

        benchWWL24(input, param, cache, accessPattern, cfg, true);
        benchPipe(input, param, cache, accessPattern, cfg, true);
        benchGinx(input, param, cache, accessPattern, cfg, true);
    }

    return 0;
}
