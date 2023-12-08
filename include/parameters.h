//
// Created by Xintong Song on 2023/12/7.
//
#ifndef HLS_YATFHE_PARAMETERS_H
#define HLS_YATFHE_PARAMETERS_H

// Parameters
// Note: BR Unfolding requires n to be divisible by the unfolding value.
#ifdef TORUS32
// LWE params
const int n = 632;
const double lwe_std_dev = 3.0517578125e-05; // 2^-15
// RLWE params
const int N = 1024, k = 1;
const double rlwe_std_dev = 2.9802322387695312e-08; // 2^-25
// RGSW params
// const int l = 6, Bg_bit = 6;
const int l = 5, Bg_bit = 5;
// KS params
const int t = 2, base_bit = 6;

#else
// From Concrete/eprint 2022/704 table 4
#define SET_4
// set 1
#if defined(SET_1)
const int n = 585, N = 1024, k = 1, Bg_bit = 8, l = 2, base_bit = 2, t = 5;
const double lwe_std_dev = 9.141776004202573E-5, rlwe_std_dev = 2.989040792967434E-8;
#elif defined(SET_2)
// set 2
const int n = 744, N = 2048, k = 1, Bg_bit = 23, l = 1, base_bit = 3 , t = 5;
const double lwe_std_dev = 7.747831515176779e-6, rlwe_std_dev = 2.2148688116005568e-16;
#elif defined(SET_3)
// set 3
const int n = 807, N = 4096, k = 1, Bg_bit = 22, l = 1, base_bit = 3, t = 5;
const double lwe_std_dev = 1.0562341599676662e-6, rlwe_std_dev = 2.168404344971009e-19;
#else
// From TFHEpp
// LWE params
const int n = 632;
const double lwe_std_dev = 3.0517578125e-05; // 2^-15
// RLWE params
const int N = 2048, k = 1;
const double rlwe_std_dev = 5.684341886080802e-14; // 2^-44
// RGSW params
// const int l = 6, Bg_bit = 6;
const int l = 4, Bg_bit = 9;
// KS params
const int t = 8, base_bit = 4;
#endif

#endif

#endif //HLS_YATFHE_PARAMETERS_H
