//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

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
    // LWE params, 132-bit
    int n {1024};
    int lweNoiseB {10}; // σ: 9.21
    int64_t qLwe{Q_32};
    int qLweBits {32};

    // RLWE params
    int torusBase {2}; // p|q
#if defined(TORUS32) // 212-bit,
    int k {1};
    int N {2048};
    int rlweNoiseB {3}; // σ: 2.21
    int64_t q {Q_32};
    int torusBits {32};
    int radixBits {4};
    int l {8};
    int lApprox {4};
#elif defined(TORUS56) // 132-bit
    int k {1};
    int N {2048};
    int rlweNoiseB {8}; // σ: 7.21
    int64_t q {Q_56};
    int torusBits {56};
    int radixBits {6};
    int l {6};
    int lApprox {3};
#else
#error "torus.h: TORUS undefined"
#endif
#endif
    int batchSize{8};
    int tasksPerThread{8};
    uint64_t qNtt {Q_60P};
    int dftBits {64};
    int driftPhase {0};

    // RGSW params
    int l2 {4}; // todo
    int lDft {dftBits / radixBits};
    // KS params: the KSK's own gadget base/level, independent of PBS's
    // radixBits/l. KS is a scalar multiply-accumulate over k*N rows, far
    // cheaper per level than a PBS external product over n iterations, so
    // it can use a much smaller base (many more, tiny digits) for a large
    // noise reduction at negligible extra cost. KSK entries are always
    // LWE_Q-domain ciphertexts (qLwe, fixed regardless of TORUS_TYPE), so
    // the usable decomposition width is capped at qLweBits.
    // Call setKsRadixBits() to change ksRadixBits.
    int ksRadixBits {2}; // b_ks
    int ksWidthBits {qLweBits}; // usable KS precision
    int ksLevel {ksWidthBits / ksRadixBits}; // ks decomposition length
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

#endif //HLS_YATFHE_YATFHE_PARAMETERS_H
