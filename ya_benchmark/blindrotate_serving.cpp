//
// Created by xsong93 on 08/25/2026.
//

#include "include/bench_out.h"
#include "include/blindrotate_serving.h"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <fstream>
#include <list>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <thread>
#include <unordered_map>

#include <sched.h>
#include <unistd.h>

#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/cache_manager.h"
#include "yautil/cache_workload_generator.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

using json = nlohmann::json;
using Clock = std::chrono::steady_clock;
using ns = std::chrono::nanoseconds;

static double toMs(const ns d) {
    return std::chrono::duration<double, std::milli>(d).count();
}

static std::string keyFile(const Method m, const int tenant) {
    switch (m) {
        case Method::Ours:
            return DiskReader::generateLazyKeyFilename(tenant);
        case Method::Tfhe:
            return DiskReader::generateGinxKeyFilename(tenant);
        default:
            return DiskReader::generateWWL24KeyFilename(tenant);
    }
}

static json summarize(std::vector<double> v) {
    json o;
    o["count"] = v.size();
    if (v.empty()) return o;
    std::sort(v.begin(), v.end());
    const double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sq = 0;
    for (double x : v) {
        sq += (x - mean) * (x - mean);
    }
    auto pct = [&](const double p) {
        auto i = static_cast<long>(std::ceil(p / 100.0 * v.size())) - 1;
        return v[std::min<long>(std::max<long>(i, 0), v.size() - 1)];
    };
    o["mean_ms"] = mean;
    o["sd_ms"] = std::sqrt(sq / v.size());
    o["p50_ms"] = pct(50);
    o["p90_ms"] = pct(90);
    o["p99_ms"] = pct(99);
    o["p99_9_ms"] = pct(99.9);
    o["max_ms"] = v.back();
    return o;
}

static std::string affinity() {
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) != 0) {
        return "unknown";
    }
    std::string out;
    int runStart = -1;
    for (int cpu = 0; cpu <= CPU_SETSIZE; ++cpu) {
        const bool in = cpu < CPU_SETSIZE && CPU_ISSET(cpu, &mask);
        if (in && runStart < 0) {
            runStart = cpu;
        }
        else if (!in && runStart >= 0) {
            if (!out.empty()) {
                out += ",";
            }
            out += std::to_string(runStart);
            if (cpu - 1 != runStart) {
                out += "-" + std::to_string(cpu - 1);
            }
            runStart = -1;
        }
    }
    return out.empty() ? "none" : out;
}

template <typename KeyT, typename LoadFn, typename RotateFn>
static json runServing(const ControlArgs& opt, const YatfheParameters& param, const std::vector<int>& pattern, const TorusPolynomial& v, const Tlwe& input, LoadFn loadAndMaybeRotate, RotateFn rotate) {
    ConcurrentKeyCache<KeyT> cache(opt.capacity);
    Scheduler queue(Scheduler::parseMode(opt.sched), opt.workers, opt.queueMax, opt.spill, opt.steal);

    std::atomic<long> shed{0};
    std::atomic<long> completed{0};
    // Device-level queue-depth instrumentation (R1 asked for SSD queue depth).
    std::atomic<long> loadsInFlight{0};
    std::atomic<long> loadsPeak{0};
    std::atomic<long long> loadSamples{0};
    std::atomic<long long> loadTickSum{0};
    std::atomic<long> devPeak{0};
    std::atomic<long long> devSamples{0};
    std::atomic<long long> devTickSum{0};
    std::vector<std::vector<Completion>> perWorker(opt.workers);

    auto makeScaled = [&] {
        auto s = std::make_unique<ScaledTlwe>(2 * param.N, param.n);
        rescaleTlweToNewMod(*s, input);
        return s;
    };

    // Warm-up
    {
        auto sTlwe = makeScaled();
        Trlwe out{param};
        for (long i = 0; i < opt.warm && i < static_cast<long>(pattern.size()); ++i) {
            const int tenant = pattern[i];
            const std::string file = keyFile(opt.method, tenant);
            bool done = false;
            auto acq = cache.acquire(tenant, [&] {
                auto k = std::make_shared<KeyT>();
                clearFileCache(file);
                done = loadAndMaybeRotate(*k, file, out, *sTlwe);
                return k;
            });
            if (!done) {
                rotate(*acq.value, out, *sTlwe);
            }
            for (int b = 1; b < opt.batch; ++b) {
                rotate(*acq.value, out, *sTlwe);
            }
        }
    }

    // Sample the kernel's in-flight read counter and the number of concurrent key loads over the measured window.
    // The first counts block-layer requests, the second counts workers loading a key.
    std::atomic<bool> monitorStop{false};
    std::thread monitor;
    const std::string inflightPath = opt.devStat.empty() ? std::string() : ("/sys/block/" + opt.devStat + "/inflight");
    if (!opt.devStat.empty()) {
        monitor = std::thread([&] {
            while (!monitorStop.load(std::memory_order_relaxed)) {
                long reads = 0, writes = 0;
                {
                    std::ifstream f(inflightPath);
                    if (f) f >> reads >> writes;
                }
                const long loading = loadsInFlight.load(std::memory_order_relaxed);
                loadTickSum.fetch_add(loading, std::memory_order_relaxed);
                loadSamples.fetch_add(1, std::memory_order_relaxed);
                devTickSum.fetch_add(reads, std::memory_order_relaxed);
                devSamples.fetch_add(1, std::memory_order_relaxed);
                long prev = devPeak.load(std::memory_order_relaxed);
                while (reads > prev && !devPeak.compare_exchange_weak(prev, reads)) {}
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        });
    }

    const auto t0 = Clock::now();
    std::thread producer([&] {
        ArrivalProcess arrivals(opt, opt.seed ^ 0x9e3779b97f4a7c15ULL);
        for (long i = 0; i < opt.requests; ++i) {
            const double offsetMs = arrivals.nextOffsetMs();
            const auto due = t0 + ns(static_cast<long long>(offsetMs * 1e6));
            std::this_thread::sleep_until(due);
            Request r{i, pattern[(opt.warm + i) % pattern.size()], Clock::now()};
            if (!queue.push(r)) {
                shed.fetch_add(1, std::memory_order_relaxed);
            }
        }
        queue.close();
    });

    std::vector<std::thread> workers;
    workers.reserve(opt.workers);
    for (int w = 0; w < opt.workers; ++w) {
        workers.emplace_back([&, w] {
            auto sTlwe = makeScaled();
            Trlwe out{param};
            auto& sink = perWorker[w];
            Request r;
            bool servedByOwner = false;
            while (queue.pop(w, r, servedByOwner)) {
                const auto picked = Clock::now();
                const std::string file = keyFile(opt.method, r.tenant);
                bool done = false;
                auto acq = cache.acquire(r.tenant, [&] {
                    const long nowLoading = loadsInFlight.fetch_add(1) + 1;
                    long prevPeak = loadsPeak.load(std::memory_order_relaxed);
                    while (nowLoading > prevPeak && !loadsPeak.compare_exchange_weak(prevPeak, nowLoading)) {}
                    auto k = std::make_shared<KeyT>();
                    clearFileCache(file);
                    done = loadAndMaybeRotate(*k, file, out, *sTlwe);
                    loadsInFlight.fetch_sub(1);
                    return k;
                });
                if (!done) {
                    rotate(*acq.value, out, *sTlwe);
                }

                // remaining bootstraps of this tenant visit reuse the loaded key:
                // this amortises the key load over a batch of size B
                for (int b = 1; b < opt.batch; ++b) {
                    rotate(*acq.value, out, *sTlwe);
                }
                const auto finished = Clock::now();
                sink.push_back(Completion {
                    r.id, r.tenant,
                    toMs(picked - r.arrival),
                    toMs(finished - picked),
                    toMs(finished - r.arrival),
                    acq.hit, acq.waitedForPeer, w, opt.batch, servedByOwner});
                completed.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }

    producer.join();
    for (auto& t : workers) {
        t.join();
    }
    const auto t1 = Clock::now();
    monitorStop.store(true, std::memory_order_relaxed);
    if (monitor.joinable()) {
        monitor.join();
    }

    std::vector<Completion> all;
    for (auto& v2 : perWorker) {
        all.insert(all.end(), v2.begin(), v2.end());
    }
    std::vector<double> sojourn, queueDelay, service, servicePerOp;
    long hits = 0, coalesced = 0, sloViolations = 0, ownedCount = 0;
    std::vector<long> perWorkerCount(opt.workers, 0);
    std::unordered_map<int, std::pair<double, long>> perTenant;
    for (const auto& c : all) {
        sojourn.push_back(c.sojournMs);
        queueDelay.push_back(c.queueMs);
        service.push_back(c.serviceMs);
        servicePerOp.push_back(c.serviceMs / std::max(1, c.ops));
        hits += c.hit;
        coalesced += c.waitedForPeer;
        sloViolations += (c.sojournMs > opt.sloMs);
        if (c.worker >= 0 && c.worker < opt.workers) {
            perWorkerCount[c.worker]++;
        }
        ownedCount += c.owned;
        auto& e = perTenant[c.tenant];
        e.first += c.sojournMs;
        e.second++;
    }

    double fair = 1.0;
    if (perTenant.size() > 1) {
        double sum = 0, sumSq = 0;
        for (auto& [t, e] : perTenant) {
            const double m = e.first / e.second;
            sum += m;
            sumSq += m * m;
        }
        fair = (sum * sum) / (perTenant.size() * sumSq);
    }
    const double wallSec =
        std::chrono::duration<double>(t1 - t0).count();

    json r;
    r["method"] = methodName(opt.method);
    r["schema_version"] = 1;
    r["config"] = {
        {"arrival", opt.arrival},
        {"offered_rate_rps", opt.rate},
        {"workers", opt.workers},
        {"requests", opt.requests},
        {"warm", opt.warm},
        {"users", opt.users},
        {"cache_capacity_keys", opt.capacity},
        {"zipf_s", opt.zipfS},
        {"seed", opt.seed},
        {"slo_ms", opt.sloMs},
        {"queue_max", opt.queueMax},
        {"batch", opt.batch},
        {"sched", opt.sched},
        {"spill_threshold", opt.spill},
        {"steal_threshold", opt.steal},
        {"burst_on_ms", opt.burstOnMs},
        {"burst_off_ms", opt.burstOffMs},
        {"burst_factor", opt.burstFactor},
        {"tag", opt.tag},
    };
    r["environment"] = {
        {"affinity_cpus", affinity()},
        {"cpus_available",
         static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN))}};
    r["throughput_rps"] = all.empty() ? 0.0 : all.size() / wallSec;
    r["offered_rps"] = opt.rate;
    r["wall_sec"] = wallSec;
    r["shed"] = shed.load();
    r["shed_rate"] = opt.requests ? static_cast<double>(shed.load()) / opt.requests : 0.0;
    r["hit_rate"] = all.empty() ? 0.0 : static_cast<double>(hits) / all.size();
    r["coalesced_loads"] = coalesced;
    r["slo_ms"] = opt.sloMs;
    r["slo_violation_rate"] = all.empty() ? 0.0 : static_cast<double>(sloViolations) / all.size();

    for (const double budget : {100.0, 300.0, 500.0}) {
        long v = 0;
        for (const auto& c : all) v += (c.sojournMs > budget);
        r["slo_violation_rate_" + std::to_string(static_cast<int>(budget))] = all.empty() ? 0.0 : static_cast<double>(v) / all.size();
    }
    r["shed_plus_completed"] = shed.load() + static_cast<long>(all.size());
    r["peak_queue_depth"] = queue.peakDepth();
    const long long nSamples = devSamples.load();
    r["load_concurrency_peak"] = loadsPeak.load();
    r["load_concurrency_mean"] = nSamples ? static_cast<double>(loadTickSum.load()) / nSamples : 0.0;
    r["device"] = opt.devStat;
    r["device_inflight_reads_peak"] = devPeak.load();
    r["device_inflight_reads_mean"] = nSamples ? static_cast<double>(devTickSum.load()) / nSamples : 0.0;
    r["device_inflight_samples"] = nSamples;
    r["batch"] = opt.batch;
    r["sched"] = opt.sched;
    r["per_worker_requests"] = perWorkerCount;
    r["tenant_fairness_jain"] = fair;
    r["locality"] = all.empty() ? 0.0 : static_cast<double>(ownedCount) / all.size();
    r["spilled"] = queue.spilled();
    r["stolen"] = queue.stolen();
    r["spill_threshold"] = opt.spill;
    r["service_per_op"] = summarize(servicePerOp);
    r["sojourn"] = summarize(sojourn);
    r["queue_delay"] = summarize(queueDelay);
    r["service"] = summarize(service);
    return r;
}

int main(int argc, char** argv) {
    ControlArgs opt;
    auto arg = [&](const char* name, const std::string& s, auto& dst) {
        const std::string pre = std::string("--") + name + "=";
        if (s.rfind(pre, 0) != 0) {
            return false;
        }
        std::istringstream(s.substr(pre.size())) >> dst;
        return true;
    };
    std::string method = "ours";
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s == "--help") {
            std::cout << "Usage: blindrotate_serving [options]\n"
                   "  --method=ours|tfhe|wwl24   scheme under test (default ours)\n"
                   "  --rate=N                   offered requests/sec (default 20)\n"
                   "  --workers=N                concurrent servers (default 4)\n"
                   "  --requests=N               measured requests (default 2000)\n"
                   "  --warm=N                   closed-loop warm-up requests (default 200)\n"
                   "  --users=N                  tenants (default 50)\n"
                   "  --cap=N                    cache slots for this method (default 10)\n"
                   "  --arrival=poisson|bursty|uniform\n"
                   "  --burston=MS --burstoff=MS --burstfactor=X\n"
                   "  --slo=MS                   SLO target for violation rate (default 100)\n"
                   "  --queuemax=N               shed beyond this depth (default 4096)\n"
                   "  --devstat=DEV              sample /sys/block/DEV/inflight (e.g. nvme0n1)\n"
                   "  --batch=B                  bootstraps per tenant visit (default 1).\n"
                   "  --sched=shared|affinity|soft\n"
                   "                             shared:   one queue, no locality, even load\n"
                   "                             affinity: hard tenant->worker hash (imbalances)\n"
                   "                             soft:     affinity + spill on admission and\n"
                   "                                       work-stealing on dequeue (recommended)\n"
                   "  --spill=N                  soft: spill when the owner lane is N deep (default 2)\n"
                   "  --steal=N                  soft: only steal from a lane N deep (default 2).\n"
                   "                             Higher preserves locality, lower favours balance.\n"
                   "  --s=F --seed=N --tag=STR\n"
                   "  --out=DIR                  directory for result json\n";
            return 0;
        }
        arg("method", s, method);
        arg("rate", s, opt.rate);
        arg("workers", s, opt.workers);
        arg("requests", s, opt.requests);
        arg("warm", s, opt.warm);
        arg("users", s, opt.users);
        arg("cap", s, opt.capacity);
        arg("arrival", s, opt.arrival);
        arg("slo", s, opt.sloMs);
        arg("queuemax", s, opt.queueMax);
        arg("s", s, opt.zipfS);
        arg("seed", s, opt.seed);
        arg("tag", s, opt.tag);
        arg("devstat", s, opt.devStat);
        {
            std::string od;
            if (arg("out", s, od))
                yabench::setBenchOutDir(od);
        }
        arg("burston", s, opt.burstOnMs);
        arg("batch", s, opt.batch);
        arg("sched", s, opt.sched);
        arg("spill", s, opt.spill);
        arg("steal", s, opt.steal);
        arg("burstoff", s, opt.burstOffMs);
        arg("burstfactor", s, opt.burstFactor);
    }
    if (method == "tfhe") {
        opt.method = Method::Tfhe;
    } else if (method == "wwl24") {
        opt.method = Method::Wwl24;
    } else {
        opt.method = Method::Ours;
    }

    YatfheParameters param{};
    initYatfhe(param);

    if (opt.batch < 1) opt.batch = 1;
    if (opt.sched != "shared" && opt.sched != "affinity" && opt.sched != "soft") {
        std::cerr << "--sched must be shared, affinity or soft\n";
        return 1;
    }
    if (opt.spill < 1) opt.spill = 1;

    std::cout << "serving: method=" << methodName(opt.method)
              << " arrival=" << opt.arrival
              << " rate=" << opt.rate << " rps"
              << " workers=" << opt.workers
              << " cap=" << opt.capacity
              << " users=" << opt.users
              << " reqs=" << opt.requests
              << " batch=" << opt.batch
              << " sched=" << opt.sched
              << " slo=" << opt.sloMs << "ms"
              << " cpus=" << affinity() << "\n";

    // Client side
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);
    vector<NttPolynomial> gdVntt;
    prepareGdV(gdVntt, v, param);

    const int p = param.torusBase;
    Integer pt = -1;
    const int slot = (pt % p + p) % p;
    Torus mu = modSwitchToTorusGeneral(slot, 2 * p, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);

    CacheWorkloadGenerator workload(opt.zipfS, opt.users, opt.seed);
    const auto pattern = workload.generateAccessPattern(opt.warm + opt.requests + 16);

    json result;
    switch (opt.method) {
        case Method::Ours:
            result = runServing<BootstrappingKeyMPLazyPipeAlt>(
                opt, param, pattern, v, input,
                [&](BootstrappingKeyMPLazyPipeAlt& k, const std::string& f, Trlwe& out, const ScaledTlwe& s) {
                    blindRotateLazyPipeAltInitNtt(out, k, s, v, f, gdVntt, param);
                    return true;
                },
                [&](const BootstrappingKeyMPLazyPipeAlt& k, Trlwe& out, const ScaledTlwe& s) {
                    blindRotateLazyPipeAltNtt(out, k, s, v, gdVntt, param);
                });
            break;
        case Method::Tfhe:
            result = runServing<BootstrappingKeyMP>(
                opt, param, pattern, v, input,
                [&](BootstrappingKeyMP& k, const std::string& f, Trlwe&, const ScaledTlwe&) {
                    deserializeBskMP(k, f, param.n);
                    return false;
                },
                [&](const BootstrappingKeyMP& k, Trlwe& out, const ScaledTlwe& s) {
                    genNoiselessTrlweSample(out, v, s);
                    blindRotateJP22Ntt(out, k, s, param);
                });
            break;
        default:
            result = runServing<BootstrappingKeyWWL24>(
                opt, param, pattern, v, input,
                [&](BootstrappingKeyWWL24& k, const std::string& f, Trlwe&, const ScaledTlwe&) {
                    deserializeBskWWL24(k, f, param.n);
                    return false;
                },
                [&](const BootstrappingKeyWWL24& k, Trlwe& out, const ScaledTlwe& s) {
                    genNoiselessTrlweSample(out, v, s);
                    blindRotateWWL24Ntt(out, k, s, param);
                });
            break;
    }

    std::string name = std::string("serving_") + (opt.method == Method::Ours ? "ours" : opt.method == Method::Tfhe ? "tfhe" : "wwl+24");
    if (!opt.tag.empty()) {
        name += "_" + opt.tag;
    }
    name += ".json";
    const std::string outPath = yabench::benchOutPath(name);
    std::ofstream out(outPath);
    out << result.dump(2) << std::endl;

    std::cout << "throughput " << result["throughput_rps"].get<double>()
              << " rps (offered " << opt.rate << ")  shed "
              << result["shed"].get<long>() << " ("
              << 100 * result["shed_rate"].get<double>() << "%)  hit "
              << result["hit_rate"].get<double>() << "  coalesced "
              << result["coalesced_loads"].get<long>() << "\n";
    std::cout << "sojourn mean " << result["sojourn"]["mean_ms"].get<double>()
              << " p99 " << result["sojourn"]["p99_ms"].get<double>()
              << " p99.9 " << result["sojourn"]["p99_9_ms"].get<double>()
              << " ms | queue mean "
              << result["queue_delay"]["mean_ms"].get<double>()
              << " p99 " << result["queue_delay"]["p99_ms"].get<double>()
              << " | service mean "
              << result["service"]["mean_ms"].get<double>() << " ms\n";
    std::cout << "service/op "
              << result["service_per_op"]["mean_ms"].get<double>()
              << " ms  fairness(Jain) "
              << result["tenant_fairness_jain"].get<double>()
              << " locality " << result["locality"].get<double>()
              << " spilled " << result["spilled"].get<long>()
              << " stolen " << result["stolen"].get<long>() << "\n";
    std::cout << "SLO " << opt.sloMs << " ms violated by "
              << 100 * result["slo_violation_rate"].get<double>() << "%"
              << "  peak queue depth "
              << result["peak_queue_depth"].get<size_t>() << "\n";
    std::cout << "wrote " << name << "\n";
    return 0;
}