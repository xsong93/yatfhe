//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_YATFHE_PARAMETERS_H
#define HLS_YATFHE_YATFHE_PARAMETERS_H

struct YatfheParameters {
// Parameters
// Note: BR Unfolding requires n to be divisible by the unfolding value.
#ifdef TORUS32
    // LWE params
    const int n = 632;
    const double lweStdDev = 3.0517578125e-05; // 2^-15
    // RLWE params
    const int N = 1024, k = 1;
    const double rlweStdDev = 2.9802322387695312e-08; // 2^-25
    // RGSW params
    // const int l = 6, Bg_bit = 6;
    const int l = 5, Bg_bit = 5;
    // KS params
    const int t = 2, baseBit = 6;

#else
// From Concrete/eprint 2022/704 table 4
#define SET_4
// set 1
#if defined(SET_1)
    const int n = 585, N = 1024, k = 1, Bg_bit = 8, l = 2, baseBit = 2, t = 5;
    const double lweStdDev = 9.141776004202573E-5, rlweStdDev = 2.989040792967434E-8;
#elif defined(SET_2)
// set 2
    const int n = 744, N = 2048, k = 1, Bg_bit = 23, l = 1, baseBit = 3 , t = 5;
    const double lweStdDev = 7.747831515176779e-6, rlweStdDev = 2.2148688116005568e-16;
#elif defined(SET_3)
// set 3
    const int n = 807, N = 4096, k = 1, Bg_bit = 22, l = 1, baseBit = 3, t = 5;
    const double lweStdDev = 1.0562341599676662e-6, rlweStdDev = 2.168404344971009e-19;
#else
    // From TFHE
    // LWE params
    int torusBits {32};
    int n {630};
    double lweStdDev {2.98023e-08}; // 2^-15
    // RLWE params
    int N {1024};
    int k {1};
    double rlweStdDev {5.684341886080802e-14}; // 2^-44
    // RGSW params
    // const int l = 6, Bg_bit = 6;
    int l {3};
    int radixBits {7}; // b
    // KS params
    int t {torusBits / radixBits}; // ks decomposition length
    int radixBase {1 << radixBits};  // 2^b
    int bHalf {radixBase / 2};
    int digitMask {radixBase - 1};
    int unfolding {1};
    int torusBase {8}; // p|q

#endif

#endif

};

#endif //HLS_YATFHE_YATFHE_PARAMETERS_H
