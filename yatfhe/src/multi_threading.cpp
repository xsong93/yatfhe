//
// Created by xintong on 4/23/25.
// Copyright (c) 2012 Jakob Progsch, Václav Zeman.
//

#include <functional>
#include "yautil/multi_threading.h"

#include <iostream>

void ThreadPool::initThreadPool() {
    auto& pool = ThreadPool::instance();
    pool.enqueue([](){ /* initialization task */ });
}

ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for(size_t i = 0; i < threads; ++i)
        workers.emplace_back([this] {
            for(;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock,
                                         [this]{ return this->stop || !this->tasks.empty(); });
                    if(this->stop && this->tasks.empty())
                        return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                task();
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