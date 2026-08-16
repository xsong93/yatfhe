//
// Created by Xintong Song on 2023/12/7.
//
#ifndef YATFHE_YATFHE_PARAMETERS_H
#define YATFHE_YATFHE_PARAMETERS_H

#include "yatfhe/torus.h"

struct YatfheParameters {
#ifdef TERNARY // ternary secret
#error "TERNARY Not Implemented"
    // // LWE params
    // int n {430};
    // int lweNoiseB {11};
    // int64_t qLwe{Q_20};
    // int qLweBits {20};
    //
    // // RLWE params
    // int k {1};
    // int N {1024};
    // int rlweNoiseB {4};
    // int64_t q {Q_27};
    // int torusBits {27};
    // int torusBase {8}; // p|q
    //
    // int radixBits {8}; // b
    // int l {4};
    // int lApprox {2};
#else // binary secret
    // LWE params, 128.3-bit
    int n {680};
    int lweNoiseB {17}; // log2 σ: 17.21
    int64_t qLwe{Q_32};
    int qLweBits {32};

    // RLWE params
    int torusBase {4}; // p|q
#if defined(TORUS32) // 212-bit,
    int k {1};
    int N {2048};
    int rlweNoiseB {3}; // log2 σ: 2.21
    int64_t q {Q_32};
    int torusBits {32};
    int radixBits {4};
    int l {8};
    int lApprox {4};
    uint64_t qNtt {Q_49P_T32};
#elif defined(TORUS56) // 132-bit
    int k {1};
    int N {2048};
    int rlweNoiseB {8}; // log2 σ: 7.21
    int64_t q {Q_56};
    int torusBits {56};
    int radixBits {14};
    int l {4};
    int lApprox {1};
    uint64_t qNtt {Q_61P_T56};
#else
#error "torus.h: TORUS undefined"
#endif
#endif
    int batchSize{8};
    int tasksPerThread{8};
    int dftBits {64};
    int driftPhase {0};

    // RGSW params
    int l2 {4}; // todo
    int lDft {dftBits / radixBits};
    int radixBase {1 << radixBits};  // 2^b
    int baseOverTwo {radixBase / 2}; // B / 2 threshold
    int digitMask {radixBase - 1};

    // KS params
    int ksRadixBits {2}; // b_ks
    int ksWidthBits {qLweBits}; // usable KS precision
    int ksLevel {ksWidthBits / ksRadixBits}; // ks decomposition length

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

    void setRadixBits(const int b) {
        radixBits = b;
        lDft = dftBits / radixBits;
        radixBase = 1 << radixBits;
        baseOverTwo = radixBase / 2;
        digitMask = radixBase - 1;
    }

    void setKsRadixBits(const int b) {
        ksRadixBits = b;
        ksWidthBits = qLweBits;
        ksLevel = ksWidthBits / ksRadixBits;
    }
};

#endif //YATFHE_YATFHE_PARAMETERS_H
