//
// Created by Xintong Song on 2023/12/7.
//
#ifndef YATFHE_TIME_COUNTER_H
#define YATFHE_TIME_COUNTER_H

#include <ctime>
#include <chrono>
#include <iostream>

using namespace std::chrono;

class TimeCounter {
private:
    time_point<steady_clock> time{};

public:
    static time_point<steady_clock> timeGlobal;

    TimeCounter();

    void resetTime();

    time_point<steady_clock> getTime();

    void printTime(const std::string& message);

    static void printTime(const std::string& message, time_point<steady_clock>& start);

    static void printTime(const std::string& message, time_point<steady_clock>& start, long interval);

    static void printTimeRefresh(const std::string& message, clock_t *start);
};

#define COUNT_TIME(MSG, CODE)                                             \
    TimeCounter::timeGlobal = steady_clock::now();               \
    CODE;                                                                 \
    TimeCounter::printTime(MSG, TimeCounter::timeGlobal);

#define BENCH_CUSTOM(MSG, CODE, COUNT)                                    \
    TimeCounter::timeGlobal = steady_clock::now();               \
    for (int benchCustom = 0; benchCustom < COUNT; benchCustom++) {       \
        CODE;                                                             \
    }                                                                     \
    TimeCounter::printTime(MSG, TimeCounter::timeGlobal, COUNT);

#define BENCH500(MSG, CODE)                                               \
    TimeCounter::timeGlobal = steady_clock::now();               \
    for (int bench500 = 0; bench500 < 500; bench500++) {                  \
        CODE;                                                             \
    }                                                                     \
    TimeCounter::printTime(MSG, TimeCounter::timeGlobal, 500);


#endif //YATFHE_TIME_COUNTER_H
