//
// Created by xsong93 on 08/25/2026.
//

#pragma once

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <list>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;

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
            if (inflight == inflight_.end()) {
                break;
            }
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
        if (index_.count(key)) {
            return;
        }
        order_.emplace_front(key, std::move(value));
        index_[key] = order_.begin();
        while (order_.size() > capacity_) {
            index_.erase(order_.back().first);
            order_.pop_back();
        }
    }

    mutable std::mutex mutex_;
    std::list<std::pair<int, std::shared_ptr<V>>> order_;
    std::unordered_map<int, typename std::list<std::pair<int, std::shared_ptr<V>>>::iterator> index_;
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
        if (s == "affinity") {
            return Mode::Affinity;
        }
        if (s == "soft") {
            return Mode::Soft;
        }
        return Mode::Shared;
    }

    Scheduler(const Mode mode, const int workers, const size_t maxDepth, const int spillThreshold, const int stealThreshold)
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
            if (mode_ == Mode::Soft && lanes_[lane].size() >= static_cast<size_t>(spill_)) {
                const int alt = shortestLaneLocked();
                if (lanes_[alt].size() + 1 < lanes_[lane].size()) {
                    lane = alt;
                    ++spilled_;
                }
            }
        }
        if (lanes_[lane].size() >= maxDepth_) {
            return false;
        }
        lanes_[lane].push_back(r);
        peak_ = std::max(peak_, lanes_[lane].size());
        cv_.notify_all();
        return true;
    }

    bool pop(const int workerId, Request& out, bool& servedByOwner) {
        std::unique_lock<std::mutex> lk(m_);
        const int own = mode_ == Mode::Shared ? 0 : workerId;
        auto ready = [&] {
            if (!lanes_[own].empty()) {
                return true;
            }
            if (mode_ == Mode::Affinity) {
                return false;
            }
            if (mode_ == Mode::Shared) {
                return anyPendingLocked();
            }
            return longestDepthLocked() >= static_cast<size_t>(steal_);
        };
        cv_.wait(lk, [&] { return closed_ || ready(); });
        const bool drainOnClose = closed_ && anyPendingLocked() && mode_ != Mode::Affinity;
        if (!ready() && !drainOnClose) {
            return false;
        }

        int lane = -1;
        if (!lanes_[own].empty()) {
            lane = own;
        } else {
            lane = longestLaneLocked();
            if (lane < 0) {
                return false;
            }
            if (mode_ == Mode::Soft) {
                ++stolen_;
            }
        }
        out = lanes_[lane].front();
        lanes_[lane].pop_front();
        servedByOwner = (mode_ != Mode::Shared) && (owner(out.tenant) == workerId);
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> g(m_);
        closed_ = true;
        cv_.notify_all();
    }

    size_t longestDepthLocked() const {
        size_t d = 0;
        for (const auto& l : lanes_) {
            d = std::max(d, l.size());
        }
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
        for (const auto& l : lanes_) {
            if (!l.empty()) {
                return true;
            }
        }
        return false;
    }
    int shortestLaneLocked() const {
        int best = 0;
        for (size_t i = 1; i < lanes_.size(); ++i) {
            if (lanes_[i].size() < lanes_[best].size()) best = static_cast<int>(i);
        }
        return best;
    }
    int longestLaneLocked() const {
        int best = -1;
        for (size_t i = 0; i < lanes_.size(); ++i) {
            if (!lanes_[i].empty() && (best < 0 || lanes_[i].size() > lanes_[best].size())) {
                best = static_cast<int>(i);
            }
        }
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

static const char* methodName(const Method m) {
    switch (m) {
        case Method::Ours: return "OURS";
        case Method::Tfhe: return "TFHE";
        default: return "WWL+24";
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
    std::string devStat;  // Block device whose kernel in-flight read counter is sampled during the run
};

class ArrivalProcess {
public:
    ArrivalProcess(const ControlArgs& o, const uint64_t seed) : opt_(o), gen_(seed), exp_(1.0) {}

    double nextOffsetMs() {
        if (opt_.arrival == "uniform") {
            cursorMs_ += 1000.0 / opt_.rate;
            return cursorMs_;
        }
        if (opt_.arrival == "bursty") {
            // Alternating phases in TIME: the on-phase lasts burstOnMs at burstFactor x rate,
            // the off-phase burstOffMs at a rate that makes the time-average equal the nominal rate.
            const double period = opt_.burstOnMs + opt_.burstOffMs;
            const double onShare = opt_.burstOnMs / period;
            const double offFactor = (1.0 - onShare * opt_.burstFactor) / std::max(1e-9, 1.0 - onShare);
            for (;;) {
                while (cursorMs_ >= phaseEndMs_) {
                    inOnPhase_ = !inOnPhase_;
                    phaseEndMs_ += inOnPhase_ ? opt_.burstOnMs : opt_.burstOffMs;
                }
                const double factor = inOnPhase_ ? opt_.burstFactor : std::max(0.05, offFactor);
                const double step = exp_(gen_) * 1000.0 / (opt_.rate * factor);
                if (cursorMs_ + step > phaseEndMs_) {
                    cursorMs_ = phaseEndMs_;  // phase has no room left: advance it
                    continue;
                }
                cursorMs_ += step;
                return cursorMs_;
            }
        }
        cursorMs_ += exp_(gen_) * 1000.0 / opt_.rate;  // Poisson
        return cursorMs_;
    }

private:
    const ControlArgs& opt_;
    std::mt19937_64 gen_;
    std::exponential_distribution<double> exp_;
    double cursorMs_{0.0};
    double phaseEndMs_{0.0};
    bool inOnPhase_{false};
};
