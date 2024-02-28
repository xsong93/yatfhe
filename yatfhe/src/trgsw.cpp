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
void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
    const auto N = param.N;
    const auto k = param.k;
    const auto l = param.l;
    const auto sigma = param.lweStdDev;
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            Trlwe& trlweSample = trgsw.trlweSamples[lvl][row];
            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[lvl][row];
            initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, N, sigma); // init b = 0 + e
            applyNtt(trlweDftSample.b, trlweSample.b);
            for (auto col = 0; col < k; col++) {
                initCoeffsViaUniformDistribution(trlweSample.a[col].coeffs, N); // init a
                applyNtt(trlweDftSample.a[col], trlweSample.a[col]);
                applyNtt(trgswKey.trlweKey.sDft[col], trgswKey.trlweKey.s[col]);
                modularAccumulate(trlweDftSample.b.coeffs, trlweDftSample.a[col].coeffs, trgswKey.trlweKey.sDft[col].coeffs);
            }
            applyIntt(trlweSample.b, trlweDftSample.b);
        }
    }
}

// output += mu * G^T
void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, int mu, const YatfheParameters& param) {
    const auto k = param.k;
    const auto l = param.l;
    const auto g = genPowersOfBgbit(param.bgBit, l);

    // add the diagonal matrix (mu * G^T)_ijk to the output
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {

            // add to a_lii
            if (row < k) {
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += mu * g[lvl]; // coeffs[0]: add mu to the constant polynomial term
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
            trgsw.trlweSamples[lvl][row].b.coeffs[0] += mu * g[lvl];
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
//            LagrangePolynomial tmp(param.N);
//            generateLagrangePolynomialWithValueAt(tmp, mu * g[layer], 0);
//            lagrangePolynomialAccumulate(output.trlweDftSamples[layer][row].b, tmp);
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






