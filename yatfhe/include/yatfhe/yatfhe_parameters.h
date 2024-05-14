//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

struct YatfheParameters {
    // LWE params
    int torusBits {32};
    long q {1 << 31};
    int torusBase {8}; // p|q
    int dftBits {64};
    int n {630};
    double lweStdDev {2.98023e-08}; // 2^-15
    // RLWE params
    static const int N {1024};
    int driftPhase {N / torusBase};
    int k {1};
    double rlweStdDev {2.98023223876953e-8}; // 2^-25
    // RGSW params
    int l {8};
    int l2 {8};
    int lDft {l * 2};
    int radixBits {4}; // b
    // KS params
    int ksLevel {torusBits / radixBits}; // ks decomposition length
    int radixBase {1 << radixBits};  // 2^b
    int baseOverTwo {radixBase / 2}; // B / 2 threshold
    int digitMask {radixBase - 1};
    int unfolding {1};
};

#endif //HLS_YATFHE_YATFHE_PARAMETERS_H
