//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_TIME_COUNTER_H
#define HLS_YATFHE_TIME_COUNTER_H

#include <ctime>
#include <chrono>
#include <iostream>
#include <mutex>

using namespace std::chrono;

class TimeCounter {
private:
    time_point<high_resolution_clock> time{};
    std::mutex mux;

public:
    TimeCounter();

    static TimeCounter* init();

    void lock();

    void unlock();

    void resetTime();

    time_point<high_resolution_clock> getTime();

    void printTime(const std::string& message);

    static void printTime(const std::string& message, time_point<high_resolution_clock>& start);

    static void printTime(const std::string& message, clock_t& start, long interval);

    static void printTimeRefresh(const std::string& message, clock_t *start);
};

#define COUNT_TIME(MSG, T, CODE) \
    T->lock();                   \
    T->resetTime();              \
    CODE;                        \
    T->printTime(MSG);           \
    T->unlock();

#endif //HLS_YATFHE_TIME_COUNTER_H
