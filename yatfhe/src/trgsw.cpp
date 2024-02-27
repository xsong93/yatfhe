//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe_parameters.h"
#include "numeric_functions.h"
#include "ntt.h"
#include "trgsw.h"
#include "trlwe.h"

// trgsw(0)
void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
    const auto N = param.N;
    const auto k = param.k;
    const auto l = param.l;
    const auto sigma = param.lweStdDev;
    for (auto i = 0; i < k + 1; i++) {
        for (auto j = 0; j < l; j++) {
            Trlwe& trlweSample = trgsw.trlweSamples[i][j];
            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[i][j];
            initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, N, sigma); // init b = 0 + e
            for (auto m = 0; m < k; m++) {
                initCoeffsViaUniformDistribution(trlweSample.a[m].coeffs, N); // init a
                applyNtt(trlweDftSample.a[m], trlweSample.a[m]);
                applyNtt(trgswKey.trlweKey.sDft[m], trgswKey.trlweKey.s[m]);
                modularMult(trlweDftSample.b.coeffs, trlweDftSample.a[m].coeffs, trgswKey.trlweKey.sDft[m].coeffs);
//                applyIntt(, trlweDftSample.b.coeffs);
            }
        }
    }
}

// trgsw(0)
void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
    const auto N = param.N;
    const auto k = param.k;
    const auto l = param.l;
    const auto sigma = param.lweStdDev;
    for (auto i = 0; i < k + 1; i++) {
        for (auto j = 0; j < l; j++) {
            Trlwe& trlweSample = trgsw.trlweSamples[i][j];
            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[i][j];
            initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, N, sigma); // init b = 0 + e
            applyNtt(trlweDftSample.b, trlweSample.b);
            for (auto m = 0; m < k; m++) {
                initCoeffsViaUniformDistribution(trlweSample.a[m].coeffs, N); // init a
                applyNtt(trlweDftSample.a[m], trlweSample.a[m]);
                applyNtt(trgswKey.trlweKey.sDft[m], trgswKey.trlweKey.s[m]);
                modularAccumulate(trlweDftSample.b.coeffs, trlweDftSample.a[m].coeffs, trgswKey.trlweKey.sDft[m].coeffs);
            }
        }
    }
}

// output += mu * H
void trgswAddInteger(Trgsw& output, int mu, const YatfheParameters& param) {
    const auto k = param.k;
    const auto l = param.l;
    const auto g = genPowersOfBgbit(param.bgBit, l);
    for (auto i = 0; i < k + 1; i++) {
        for (auto j = 0; j < l; j++) {
            output.trlweSamples[i][j].a[i].coeffs[0] += mu * g[j]; // coeffs[0]: add mu to the constant term
        }
    }
}

//void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const YatfheParameters& param) {
//    trgswKey.trlweKey = trlweKey;
//    trgswKey.l = param.l;
//    trgswKey.bgBit = param.bgBit;
//}
//
//void initTrgswSample(Trgsw& trgsw, const YatfheParameters& param) {
//    const int l = param.l;
//    const int k = param.k;
//    const int N = param.N;
//    const int kpl = l * (k + 1);
////    trgsw.trlweSamples.resize(kpl);
//    for (int i = 0; i < kpl; i++) {
//        initTrlweSample(trgsw.trlweSamples[i], k, N);
//    }
//}
//
//void initTrgswDftSample(TrgswDft& trgswDftSample, const YatfheParameters& param) {
//    const int l = param.l;
//    const int k = param.k;
//    const int N = param.N;
//    const int kpl = l *  (k + 1);
////    trgswDftSample.trlweDftSamples.resize(kpl);
//    for (int i = 0; i < kpl; i++) {
//        initTrlweDftSample(trgswDftSample.trlweDftSamples[i], k, N);
//    }
//}

//void genNoiselessTrgswSample(Trgsw &trgswSample, Torus msg, const YatfheParameters& param) {
//    const int l = param.l;
//    const int bgBit = param.bgBit;
//    const int k = param.k;
//    const int N = param.N;
//    const int kpl = (k + 1) * l;
//    for (int i = 0; i < l; i++) {
//        const uint64_t h = 1UL << (sizeof(Torus) * 8 - (i + 1) * bgBit);
//        for (int j = 0; j < k; j++) {
//            trgswSample.trlweSamples[j * l + i].a[j].coeffs[0] += msg * h;
//        }
//        trgswSample.trlweSamples[k * l + i].b.coeffs[0] += msg * h;
//    }
//}

//void deleteTrgswKey(TrgswKey& trgswKey) {
////    trgswKey.trlweKey = nullptr;
//}






