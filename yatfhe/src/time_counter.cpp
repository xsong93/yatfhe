//
// Created by Xintong Song on 2023/12/8.
//

#include <iostream>
#include "time_counter.h"
#include "control_helper.h"

TimeCounter::TimeCounter() {
    time = high_resolution_clock::now();
    printTime("init timer");
}

void TimeCounter::lock () {
    mux.lock();
}

void TimeCounter::unlock () {
    mux.unlock();
}

void TimeCounter::resetTime() {
    time = high_resolution_clock::now();
}

time_point<high_resolution_clock> TimeCounter::getTime() {
    return time;
}

void TimeCounter::printTime(const std::string& message) {
#ifdef ENABLE_TIMER
    std::cout << "elapsed time (" << message << ") in us: " << duration_cast<microseconds>(high_resolution_clock::now() - getTime()).count() << std::endl;
#endif
}

void TimeCounter::printTime(const std::string& message, time_point<high_resolution_clock>& start) {
#ifdef ENABLE_TIMER
    std::cout << "elapsed time (" << message << ") in us: " << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << std::endl;
#endif
}

void TimeCounter::printTime(const std::string& message, clock_t& start, long interval) {
#ifdef ENABLE_TIMER
    if (interval <= 0) {
        interval = 1;
    }
    std::cout << "elapsed time (" << message << ") in us: " << (clock() - start) / interval << std::endl;
#endif
}

void TimeCounter::printTimeRefresh(const std::string& message, clock_t *start) {
#ifdef ENABLE_TIMER
    std::cout << "elapsed time (" << message << ") in us: " << (clock() - *start) << std::endl;
    *start = clock();
#endif
}

