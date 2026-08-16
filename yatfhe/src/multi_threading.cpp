//
// Created by xintong on 4/23/25.
// Copyright (c) 2012 Jakob Progsch, Václav Zeman.
//

#include <functional>
#include "yautil/multi_threading.h"

#include <iostream>

#ifdef __linux__
#include <sched.h>
#endif

namespace {
    thread_local bool insideWorker = false;
}

bool ThreadPool::onWorkerThread() { return insideWorker; }

unsigned ThreadPool::usableConcurrency() {
#ifdef __linux__
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int allowed = CPU_COUNT(&set);
        if (allowed > 0) return static_cast<unsigned>(allowed);
    }
#endif
    const unsigned reported = std::thread::hardware_concurrency();
    return reported > 0 ? reported : 1; // 0 means unknown, and would hang the pool
}

void ThreadPool::initThreadPool() {
    auto& pool = ThreadPool::instance();
    pool.enqueue([](){ /* initialization task */ });
}

void TaskGroup::finish(std::exception_ptr error) {
    if (error) {
        std::lock_guard<std::mutex> lock(mutex);
        if (!firstError) firstError = error;
    }
    if (remaining.fetch_sub(1, std::memory_order_acq_rel) == 1) {
        std::lock_guard<std::mutex> lock(mutex);
        condition.notify_all();
    }
}

void TaskGroup::wait() {
    std::unique_lock<std::mutex> lock(mutex);
    condition.wait(lock, [this] { return remaining.load(std::memory_order_acquire) == 0; });
    if (firstError) {
        const auto error = firstError;
        firstError = nullptr;
        std::rethrow_exception(error);
    }
}

ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for(size_t i = 0; i < threads; ++i)
        workers.emplace_back([this] {
            insideWorker = true;
            for(;;) {
                Job job;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock,
                                         [this]{ return this->stop || !this->tasks.empty(); });
                    if(this->stop && this->tasks.empty())
                        return;
                    job = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                if (job.indexed) {
                    std::exception_ptr error;
                    try {
                        job.indexed(job.callable, job.index);
                    } catch (...) {
                        error = std::current_exception();
                    }
                    job.group->finish(error);   // must run even on throw, or wait() hangs
                } else {
                    job.standalone();
                }
            }
        });
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
    std::cout << "thead end" << std::endl;
}