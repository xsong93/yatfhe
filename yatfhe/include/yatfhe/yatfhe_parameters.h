//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

#include "yatfhe/torus.h"
#include "yautil/control_helper.h"

struct YatfheParameters {
#ifdef TERNARY // ternary secret, 128-bit
    // LWE params
    int n {430};
    double lweStdDev {9.5367431640625e-07}; // 2^-20
    int64_t qLwe{Q_20};

    // RLWE params
    int k {1};
    int N {1024};
    double rlweStdDev {2.9802322387695312e-08}; // 2^-25
    int64_t q {Q_27};
    int torusBits {27};
    int torusBase {8}; // p|q

    int radixBits {8}; // b
    int l {4};
    int lApprox {2};
#else // binary secret, 128-bit
    // LWE params
    int n {512};
//    double lweStdDev {9.5367431640625e-07}; // 2^-20
    double lweStdDev {0.001953125}; // 2^-9
    int64_t qLwe{Q_20};

    // RLWE params, 128-bit
    int k {1};
    int N {2048};
    // double rlweStdDev {7.275957614183426e-12}; // 2^-37
    double rlweStdDev {4.656612873077393e-10}; // 2^-31
    int64_t q {Q_32};
    int torusBits {32};
    int torusBase {8}; // p|q

    int radixBits {8}; // b
    int l {4};
    int lApprox {2};
#endif
    int batchSize{8};
    int tasksPerThread{8};
    uint64_t qNtt {Q_50P};
    int dftBits {64};
    int driftPhase {2*N/n};

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
