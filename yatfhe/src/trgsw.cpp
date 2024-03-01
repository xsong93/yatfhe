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






