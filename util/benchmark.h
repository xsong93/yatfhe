//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_BENCHMARK_H
#define HLS_YATFHE_BENCHMARK_H
#include <sys/time.h>
#include <iostream>

uint64_t get_time() {
    struct timeval tv{};
    gettimeofday(&tv, nullptr);
    return (tv.tv_usec) + (tv.tv_sec * 1000000);
}

#define MAX_EXECS 10000
uint64_t __g_clock_begin, __g_clock_end, __g_clock_array[MAX_EXECS];

void print_bench(char * msg, int itr) {
    uint64_t mean = 0;
    uint64_t sq_err = 0;
    for (int i = 0; i < itr; i++) {
        mean += __g_clock_array[i];
    }
    mean /= itr;
    for (int i = 0; i < itr; i++) {
        sq_err += (__g_clock_array[i] - mean) * (__g_clock_array[i] - mean);
    }
    double stddev = sqrt(((double) sq_err) / itr);
    printf("%s: |%llu,%03llu,%03llu|μs +- |%lf\n", msg, mean/1000000, (mean/1000)%1000, mean%1000, stddev);
}

#define BENCHMARK(NAME, REP, MSG, CODE) \
    for (size_t ___i = 0; ___i < REP; ___i++) { \
        __g_clock_begin = get_time(); \
        CODE; \
        __g_clock_end = get_time(); \
        __g_clock_array[___i] = __g_clock_end - __g_clock_begin; \
    }\
    print_bench(MSG, REP);

#endif //HLS_YATFHE_BENCHMARK_H
