//
// Created for anonymous review.
//

#include <iostream>
#include "yautil/time_counter.h"
#include "yautil/control_helper.h"
#include "yautil/tool.h"

time_point<steady_clock> TimeCounter::timeGlobal = steady_clock::now();

TimeCounter::TimeCounter() {
    time = steady_clock::now();
    printTime("init timer");
}

void TimeCounter::resetTime() {
    time = steady_clock::now();
}

time_point<steady_clock> TimeCounter::getTime() {
    return time;
}

void TimeCounter::printTime(const std::string& message) {
#ifdef ENABLE_TIMER
    printMsg(duration_cast<microseconds>(steady_clock::now() - getTime()).count(),
             "elapsed time (" + message + ") in us");
#endif
}

void TimeCounter::printTime(const std::string& message, time_point<steady_clock>& start) {
#ifdef ENABLE_TIMER
    printMsg(duration_cast<microseconds>(steady_clock::now() - start).count(),
             "elapsed time (" + message + ") in us");
#endif
}

void TimeCounter::printTime(const std::string& message, time_point<steady_clock>& start, long interval) {
#ifdef ENABLE_TIMER
    if (interval <= 0) {
        interval = 1;
    }
    printMsg((duration_cast<microseconds>(steady_clock::now() - start).count()) / interval,
             "elapsed time (" + message + ") in us");
#endif
}

void TimeCounter::printTimeRefresh(const std::string& message, clock_t *start) {
#ifdef ENABLE_TIMER
    std::cout << "elapsed time (" << message << ") in us: " << (clock() - *start) << std::endl;
    *start = clock();
#endif
}

