//
// Created by xintong on 4/23/25.
// Copyright (c) 2012 Jakob Progsch, Václav Zeman.
//

#ifndef HLS_YATFHE_MULTI_THREADING_H
#define HLS_YATFHE_MULTI_THREADING_H

#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>

class ThreadPool {
public:

    static void initThreadPool();

    static ThreadPool& instance() {
        static ThreadPool instance(std::thread::hardware_concurrency() / 2);
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
            tasks.emplace([task](){ (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ~ThreadPool();

private:
    explicit ThreadPool(size_t threads);
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

#endif //HLS_YATFHE_MULTI_THREADING_H
