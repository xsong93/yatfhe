//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

#include "yatfhe/torus.h"
#include "yautil/control_helper.h"

struct YatfheParameters {
#ifdef PAR_SET1 // 128-bit
    int n {586};
    int k {2};
    int N {512};
    int radixBits {8}; // b
    int l {2};
#else //110-bit
    int n {500};
    int k {1};
    int N {1024};
    int radixBits {10}; // b
    int l {2};
#endif
    // LWE params
    int torusBits {32};
    int64_t q {Q_32};
    int torusBase {8}; // p|q
    int dftBits {64};
    double lweStdDev {3.0517578125e-05}; // 2^-15
    // RLWE params
    int driftPhase {N / torusBase / 2};
    double rlweStdDev {2.9802322387695312e-08}; // 2^-25
    // RGSW params
    int l2 {4}; // todo
    int lDft {dftBits / radixBits};
    // KS params
    int ksLevel {torusBits / radixBits}; // ks decomposition length
    int radixBase {1 << radixBits};  // 2^b
    int baseOverTwo {radixBase / 2}; // B / 2 threshold
    int digitMask {radixBase - 1};
    int unfolding {1};
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
