//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

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
    long q {1 << 31};
    int torusBase {8}; // p|q
    int dftBits {64};
    double lweStdDev {2.98023e-08}; // 2^-15
    // RLWE params
    int driftPhase {N / torusBase};
    double rlweStdDev {2.98023223876953e-8}; // 2^-25
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
    int d {4}; // # of primes
    int dh {2};
    int dl{2};
    int qd[4] {251, 241, 239, 233};
    int qh[2] {qd[0], qd[1]};
    int ql[2] {qd[2], qd[3]};
    int qHigh {qd[0] * qd[1]};
    int qLow {qd[2] * qd[3]};
    long qCRT{qd[0] * qd[1] * qd[2] * qd[3]};

};

#endif //HLS_YATFHE_YATFHE_PARAMETERS_H
