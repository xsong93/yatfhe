//
// Created by xsong93 on 08/25/2026.
//

#include "include/bench_out.h"
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

namespace {

// Shared key cache
template <typename V>
class ConcurrentKeyCache {
public:
    explicit ConcurrentKeyCache(const size_t capacity) : capacity_(capacity) {}

    struct Acquired {
        std::shared_ptr<V> value;
        bool hit{};
        bool loadedHere{};
        bool waitedForPeer{};
    };

    template <typename Loader>
    Acquired acquire(const int key, Loader&& loader) {
        std::unique_lock<std::mutex> lock(mutex_);
        for (;;) {
            if (const auto it = index_.find(key); it != index_.end()) {
                order_.splice(order_.begin(), order_, it->second);
                return {it->second->second, true, false, false};
            }
            const auto inflight = inflight_.find(key);
            if (inflight == inflight_.end()) break;
            auto state = inflight->second;
            ++state->waiters;
            state->cv.wait(lock, [&] { return state->done; });
            if (state->value) {
                return {state->value, false, false, true};
            }
        }

        auto state = std::make_shared<LoadState>();
        inflight_.emplace(key, state);
        lock.unlock();

        std::shared_ptr<V> loaded;
        try {
            loaded = loader();
        } catch (...) {
            lock.lock();
            state->done = true;
            inflight_.erase(key);
            state->cv.notify_all();
            throw;
        }

        lock.lock();
        insertLocked(key, loaded);
        state->value = loaded;
        state->done = true;
        inflight_.erase(key);
        state->cv.notify_all();
        return {loaded, false, true, false};
    }

    size_t size() const {
        std::lock_guard<std::mutex> g(mutex_);
        return order_.size();
    }

private:
    struct LoadState {
        std::condition_variable cv;
        std::shared_ptr<V> value;
        bool done{false};
        int waiters{0};
    };

    void insertLocked(const int key, std::shared_ptr<V> value) {
        if (index_.count(key)) return;
        order_.emplace_front(key, std::move(value));
        index_[key] = order_.begin();
        while (order_.size() > capacity_) {
            index_.erase(order_.back().first);
            order_.pop_back();
        }
    }

    mutable std::mutex mutex_;
    std::list<std::pair<int, std::shared_ptr<V>>> order_;
    std::unordered_map<int,
                       typename std::list<std::pair<int, std::shared_ptr<V>>>::iterator>
        index_;
    std::unordered_map<int, std::shared_ptr<LoadState>> inflight_;
    size_t capacity_;
};

struct Request {
    long id{};
    int tenant{};
    Clock::time_point arrival;
};

struct Completion {
    long id{};
    int tenant{};
    double queueMs{};
    double serviceMs{};
    double sojournMs{};
    bool hit{};
    bool waitedForPeer{};
    int worker{};
    int ops{};
    bool owned{};
};

class Scheduler {
public:
    enum class Mode { Shared, Affinity, Soft };

    static Mode parseMode(const std::string& s) {
        if (s == "affinity") return Mode::Affinity;
        if (s == "soft") return Mode::Soft;
        return Mode::Shared;
    }

    Scheduler(const Mode mode, const int workers, const size_t maxDepth,
              const int spillThreshold, const int stealThreshold)
        : mode_(mode),
          workers_(workers),
          maxDepth_(maxDepth),
          spill_(spillThreshold),
          steal_(stealThreshold),
          lanes_(mode == Mode::Shared ? 1 : workers) {}

    bool push(const Request& r) {
        std::lock_guard<std::mutex> g(m_);
        int lane = 0;
        if (mode_ != Mode::Shared) {
            lane = owner(r.tenant);
            if (mode_ == Mode::Soft && lanes_[lane].size() >=
                static_cast<size_t>(spill_)) {
                const int alt = shortestLaneLocked();
                if (lanes_[alt].size() + 1 < lanes_[lane].size()) {
                    lane = alt;
                    ++spilled_;
                }
            }
        }
        if (lanes_[lane].size() >= maxDepth_) return false;
        lanes_[lane].push_back(r);
        peak_ = std::max(peak_, lanes_[lane].size());
        cv_.notify_all();
        return true;
    }

    bool pop(const int workerId, Request& out, bool& servedByOwner) {
        std::unique_lock<std::mutex> lk(m_);
        const int own = mode_ == Mode::Shared ? 0 : workerId;
        auto ready = [&] {
            if (!lanes_[own].empty()) return true;
            if (mode_ == Mode::Affinity) return false;
            if (mode_ == Mode::Shared) return anyPendingLocked();
            return longestDepthLocked() >= static_cast<size_t>(steal_);
        };
        cv_.wait(lk, [&] { return closed_ || ready(); });
        if (!ready()) {
            if (closed_ && anyPendingLocked() && mode_ != Mode::Affinity) {
            } else {
                return false;
            }
        }

        int lane = -1;
        if (!lanes_[own].empty()) {
            lane = own;
        } else {
            lane = longestLaneLocked();
            if (lane < 0) return false;
            if (mode_ == Mode::Soft) ++stolen_;
        }
        out = lanes_[lane].front();
        lanes_[lane].pop_front();
        servedByOwner = (mode_ != Mode::Shared) &&
                        (owner(out.tenant) == workerId);
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> g(m_);
        closed_ = true;
        cv_.notify_all();
    }

    size_t longestDepthLocked() const {
        size_t d = 0;
        for (const auto& l : lanes_) d = std::max(d, l.size());
        return d;
    }
    size_t peakDepth() const { std::lock_guard<std::mutex> g(m_); return peak_; }
    long spilled() const { std::lock_guard<std::mutex> g(m_); return spilled_; }
    long stolen() const { std::lock_guard<std::mutex> g(m_); return stolen_; }

private:
    int owner(const int tenant) const {
        return static_cast<int>(static_cast<unsigned>(tenant) % workers_);
    }
    bool anyPendingLocked() const {
        for (const auto& l : lanes_) if (!l.empty()) return true;
        return false;
    }
    int shortestLaneLocked() const {
        int best = 0;
        for (size_t i = 1; i < lanes_.size(); ++i)
            if (lanes_[i].size() < lanes_[best].size()) best = static_cast<int>(i);
        return best;
    }
    int longestLaneLocked() const {
        int best = -1;
        for (size_t i = 0; i < lanes_.size(); ++i)
            if (!lanes_[i].empty() &&
                (best < 0 || lanes_[i].size() > lanes_[best].size()))
                best = static_cast<int>(i);
        return best;
    }

    mutable std::mutex m_;
    std::condition_variable cv_;
    Mode mode_;
    int workers_;
    size_t maxDepth_;
    int spill_;
    int steal_;
    std::vector<std::deque<Request>> lanes_;
    size_t peak_{0};
    long spilled_{0};
    long stolen_{0};
    bool closed_{false};
};

enum class Method { Ours, Tfhe, Wwl24 };

}  // namespace

static const char* methodName(const Method m) {
    switch (m) {
        case Method::Ours: return "OURS";
        case Method::Tfhe: return "TFHE";
        default: return "WWL+24";
    }
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

struct ControlArgs {
    Method method{Method::Ours};
    double rate{20.0};
    int workers{4};
    long requests{2000};
    int users{50};
    int capacity{10};
    double zipfS{0.83};
    uint64_t seed{1};
    std::string arrival{"poisson"};
    double burstOnMs{200.0};
    double burstOffMs{800.0};
    double burstFactor{4.0};
    double sloMs{100.0};
    size_t queueMax{4096};
    long warm{200};
    int batch{1};
    int spill{2};
    int steal{2};
    std::string sched{"shared"};
    std::string tag;
};

class ArrivalProcess {
public:
    ArrivalProcess(const ControlArgs& o, const uint64_t seed)
        : opt_(o), gen_(seed), exp_(1.0) {}

    double nextOffsetMs() {
        if (opt_.arrival == "uniform") {
            cursorMs_ += 1000.0 / opt_.rate;
            return cursorMs_;
        }
        if (opt_.arrival == "bursty") {
            advanceBurstPhase();
            const double period = opt_.burstOnMs + opt_.burstOffMs;
            const double onShare = opt_.burstOnMs / period;
            const double offFactor =
                (1.0 - onShare * opt_.burstFactor) / std::max(1e-9, 1.0 - onShare);
            const double factor =
                inOnPhase_ ? opt_.burstFactor : std::max(0.05, offFactor);
            cursorMs_ += exp_(gen_) * 1000.0 / (opt_.rate * factor);
            return cursorMs_;
        }
        cursorMs_ += exp_(gen_) * 1000.0 / opt_.rate;  // Poisson
        return cursorMs_;
    }

private:
    void advanceBurstPhase() {
        const double period = opt_.burstOnMs + opt_.burstOffMs;
        const double phase = std::fmod(cursorMs_, period);
        inOnPhase_ = phase < opt_.burstOnMs;
    }

    const ControlArgs& opt_;
    std::mt19937_64 gen_;
    std::exponential_distribution<double> exp_;
    double cursorMs_{0.0};
    bool inOnPhase_{true};
};

static json summarize(std::vector<double> v) {
    json o;
    o["count"] = v.size();
    if (v.empty()) return o;
    std::sort(v.begin(), v.end());
    const double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sq = 0;
    for (double x : v) sq += (x - mean) * (x - mean);
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
    if (sched_getaffinity(0, sizeof(mask), &mask) != 0) return "unknown";
    std::string out;
    int runStart = -1;
    for (int cpu = 0; cpu <= CPU_SETSIZE; ++cpu) {
        const bool in = cpu < CPU_SETSIZE && CPU_ISSET(cpu, &mask);
        if (in && runStart < 0) runStart = cpu;
        else if (!in && runStart >= 0) {
            if (!out.empty()) out += ",";
            out += std::to_string(runStart);
            if (cpu - 1 != runStart) out += "-" + std::to_string(cpu - 1);
            runStart = -1;
        }
    }
    return out.empty() ? "none" : out;
}

template <typename KeyT, typename LoadFn, typename RotateFn>
static json runServing(const ControlArgs& opt, const YatfheParameters& param,
                       const std::vector<int>& pattern,
                       const TorusPolynomial& v, const Tlwe& input,
                       LoadFn loadAndMaybeRotate, RotateFn rotate) {
    ConcurrentKeyCache<KeyT> cache(opt.capacity);
    Scheduler queue(Scheduler::parseMode(opt.sched), opt.workers, opt.queueMax,
                    opt.spill, opt.steal);

    std::atomic<long> shed{0};
    std::atomic<long> completed{0};
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
        for (long i = 0; i < opt.warm && i < static_cast<long>(pattern.size());
             ++i) {
            const int tenant = pattern[i];
            const std::string file = keyFile(opt.method, tenant);
            bool done = false;
            auto acq = cache.acquire(tenant, [&] {
                auto k = std::make_shared<KeyT>();
                clearFileCache(file);
                done = loadAndMaybeRotate(*k, file, out, *sTlwe);
                return k;
            });
            if (!done) rotate(*acq.value, out, *sTlwe);
            for (int b = 1; b < opt.batch; ++b) rotate(*acq.value, out, *sTlwe);
        }
    }

    const auto t0 = Clock::now();
    std::thread producer([&] {
        ArrivalProcess arrivals(opt, opt.seed ^ 0x9e3779b97f4a7c15ULL);
        for (long i = 0; i < opt.requests; ++i) {
            const double offsetMs = arrivals.nextOffsetMs();
            const auto due = t0 + ns(static_cast<long long>(offsetMs * 1e6));
            std::this_thread::sleep_until(due);
            Request r{i, pattern[(opt.warm + i) % pattern.size()], Clock::now()};
            if (!queue.push(r)) shed.fetch_add(1, std::memory_order_relaxed);
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
                    auto k = std::make_shared<KeyT>();
                    clearFileCache(file);
                    done = loadAndMaybeRotate(*k, file, out, *sTlwe);
                    return k;
                });
                if (!done) rotate(*acq.value, out, *sTlwe);
                // remaining bootstraps of this tenant visit reuse the loaded key:
                // this amortises the key load over a batch of size B
                for (int b = 1; b < opt.batch; ++b)
                    rotate(*acq.value, out, *sTlwe);
                const auto finished = Clock::now();
                sink.push_back(Completion{
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
    for (auto& t : workers) t.join();
    const auto t1 = Clock::now();

    std::vector<Completion> all;
    for (auto& v2 : perWorker) all.insert(all.end(), v2.begin(), v2.end());
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
        if (c.worker >= 0 && c.worker < opt.workers)
            perWorkerCount[c.worker]++;
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
    r["shed_rate"] = opt.requests
                         ? static_cast<double>(shed.load()) / opt.requests
                         : 0.0;
    r["hit_rate"] = all.empty() ? 0.0
                                : static_cast<double>(hits) / all.size();
    r["coalesced_loads"] = coalesced;
    r["slo_ms"] = opt.sloMs;
    r["slo_violation_rate"] = all.empty()
                                  ? 0.0
                                  : static_cast<double>(sloViolations) /
                                        all.size();
    r["peak_queue_depth"] = queue.peakDepth();
    r["batch"] = opt.batch;
    r["sched"] = opt.sched;
    r["per_worker_requests"] = perWorkerCount;
    r["tenant_fairness_jain"] = fair;
    r["locality"] = all.empty() ? 0.0
                                : static_cast<double>(ownedCount) / all.size();
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
        if (s.rfind(pre, 0) != 0) return false;
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
    if (method == "tfhe")
        opt.method = Method::Tfhe;
    else if (method == "wwl24")
        opt.method = Method::Wwl24;
    else
        opt.method = Method::Ours;

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
                [&](BootstrappingKeyMPLazyPipeAlt& k, const std::string& f,
                    Trlwe& out, const ScaledTlwe& s) {
                    blindRotateLazyPipeAltInitNtt(out, k, s, v, f, param);
                    return true;
                },
                [&](const BootstrappingKeyMPLazyPipeAlt& k, Trlwe& out,
                    const ScaledTlwe& s) {
                    blindRotateLazyPipeAltNtt(out, k, s, v, param);
                });
            break;
        case Method::Tfhe:
            result = runServing<BootstrappingKeyMP>(
                opt, param, pattern, v, input,
                [&](BootstrappingKeyMP& k, const std::string& f, Trlwe&,
                    const ScaledTlwe&) {
                    deserializeBskMP(k, f, param.n);
                    return false;
                },
                [&](const BootstrappingKeyMP& k, Trlwe& out,
                    const ScaledTlwe& s) {
                    genNoiselessTrlweSample(out, v, s);
                    blindRotateJP22Ntt(out, k, s, param);
                });
            break;
        default:
            result = runServing<BootstrappingKeyWWL24>(
                opt, param, pattern, v, input,
                [&](BootstrappingKeyWWL24& k, const std::string& f, Trlwe&,
                    const ScaledTlwe&) {
                    deserializeBskWWL24(k, f, param.n);
                    return false;
                },
                [&](const BootstrappingKeyWWL24& k, Trlwe& out,
                    const ScaledTlwe& s) {
                    genNoiselessTrlweSample(out, v, s);
                    blindRotateWWL24Ntt(out, k, s, param);
                });
            break;
    }

    std::string name = std::string("serving_") +
                       (opt.method == Method::Ours
                            ? "ours"
                            : opt.method == Method::Tfhe ? "tfhe" : "wwl+24");
    if (!opt.tag.empty()) name += "_" + opt.tag;
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