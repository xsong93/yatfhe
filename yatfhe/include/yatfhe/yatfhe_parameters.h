//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

#include "yatfhe/torus.h"

struct YatfheParameters {
#ifdef TERNARY // ternary secret, 128-bit
    // LWE params
    int n {430};
    int lweNoiseB {11};
    int64_t qLwe{Q_20};

    // RLWE params
    int k {1};
    int N {1024};
    int rlweNoiseB {4};
    int64_t q {Q_27};
    int torusBits {27};
    int torusBase {8}; // p|q

    int radixBits {8}; // b
    int l {4};
    int lApprox {2};
#else // binary secret, 128-bit
    // LWE params
    int n {512};
    int lweNoiseB {23};
    int64_t qLwe{Q_32};

    // RLWE params, 128-bit
    int k {1};
    int N {2048};
#if defined(TORUS32)
    int rlweNoiseB {3};
    int64_t q {Q_32};
    int torusBits {32};
    int radixBits {4};
    int l {8};
    int lApprox {4};
#elif defined(TORUS36)
    int rlweNoiseB {1};
    int64_t q {Q_36};
    int torusBits {36};
    int radixBits {6};
    int l {6};
    int lApprox {4};
#elif defined(TORUS40)
    int rlweNoiseB {1};
    int64_t q {Q_40};
    int torusBits {40};
    int radixBits {10};
    int l {4};
    int lApprox {3};
#elif defined(TORUS42)
    int rlweNoiseB {3};
    int64_t q {Q_42};
    int torusBits {42};
    int radixBits {6};
    int l {7};
    int lApprox {5};
#else
#error "torus.h: TORUS undefined"
#endif
#endif
    int torusBase {4}; // p|q
    int batchSize{8};
    int tasksPerThread{8};
    uint64_t qNtt {Q_50P};
    int dftBits {64};
    int driftPhase {0};

    // RGSW params
    int l2 {4}; // todo
    int lDft {dftBits / radixBits};
    // KS params
    int ksLevel {torusBits / radixBits}; // ks decomposition length
    int radixBase {1 << radixBits};  // 2^b
    int baseOverTwo {radixBase / 2}; // B / 2 threshold
    int digitMask {radixBase - 1};
    int group {1};
    // CRT params
    int d{NUM_PRIMES}; // # of primes
    int dh{NUM_HIGH_PRIMES};
    int dl{NUM_LOW_PRIMES};
    int qd[NUM_PRIMES] {QD_CRT[0], QD_CRT[1], QD_CRT[2], QD_CRT[3]};
    int qdHalf[NUM_PRIMES] {qd[0]/2, qd[1]/2, qd[2]/2, qd[3]/2};
    int qh[NUM_HIGH_PRIMES] {qd[0], qd[1]};
    int ql[NUM_LOW_PRIMES] {qd[2], qd[3]};
    int qHigh {qd[0] * qd[1]};
    int qLow {qd[2] * qd[3]};
    long qCRT{static_cast<long>(qHigh) * static_cast<long>(qLow)};
    int qLowDivQl[NUM_LOW_PRIMES] {qLow / ql[0], qLow / ql[1]};
    int taoU[NUM_PRIMES] {1, 1, 1, 1};
    int taoUInv[NUM_PRIMES] {1, 1, 1, 1};
    long w[NUM_HIGH_PRIMES] {1, 1};
    long z[NUM_PRIMES] {1, 1, 1, 1};
};

#endif //HLS_YATFHE_YATFHE_PARAMETERS_H
