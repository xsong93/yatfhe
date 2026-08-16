//
// Created by xintong on 4/23/25.
// Copyright (c) 2012 Jakob Progsch, Václav Zeman.
//

#ifndef YATFHE_MULTI_THREADING_H
#define YATFHE_MULTI_THREADING_H

#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <atomic>
#include <exception>

// Counter-based join for a batch of tasks that are fired together.
// enqueue() costs a packaged_task allocation plus a condition_variable per task.
// A group waits on one counter.
class TaskGroup {
public:
    void wait();

private:
    friend class ThreadPool;
    void expect(int count) { remaining.fetch_add(count, std::memory_order_relaxed); }
    void finish(std::exception_ptr error);

    std::atomic<int> remaining{0};
    std::mutex mutex;
    std::condition_variable condition;
    std::exception_ptr firstError;
};

class ThreadPool {
public:

    static void initThreadPool();

    // CPUs this process may actually run on.
    static unsigned usableConcurrency();

    // True while running inside a pool worker.
    static bool onWorkerThread();

    static ThreadPool& instance() {
        static ThreadPool instance(usableConcurrency());
        return instance;
    }

    static std::vector<std::future<void>>& getFutures(size_t batchSize) {
        static std::vector<std::future<void>> futures;
        futures.reserve(batchSize);
        return futures;
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<std::result_of_t<F(Args...)>> {
        using return_type = std::result_of_t<F(Args...)>;
        auto task = std::make_shared<std::packaged_task<return_type()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");
            tasks.push(Job{[task](){ (*task)(); }});
        }
        condition.notify_one();
        return res;
    }

    // Run fn(0), ..., fn(count-1) on the pool, joining through group.wait().
    template<class F>
    void run(TaskGroup& group, const int count, const F& fn) {
        if (count <= 0) return;
        group.expect(count);
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if(stop)
                throw std::runtime_error("run on stopped ThreadPool");
            for (int i = 0; i < count; ++i)
                tasks.push(Job{{}, &invokeIndexed<F>, &fn, i, &group});
        }
        for (int i = 0; i < count; ++i)
            condition.notify_one();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ~ThreadPool();

private:
    struct Job {
        std::function<void()> standalone;
        void (*indexed)(const void*, int) {nullptr};
        const void* callable {nullptr};
        int index {0};
        TaskGroup* group {nullptr};
    };

    template<class F>
    static void invokeIndexed(const void* callable, const int index) {
        (*static_cast<const F*>(callable))(index);
    }

    explicit ThreadPool(size_t threads);
    std::vector<std::thread> workers;
    std::queue<Job> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

#endif //YATFHE_MULTI_THREADING_H
