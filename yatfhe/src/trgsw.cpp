//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"
#include "yatfhe/ntt.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"

void trgswEncrypt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey, const Integer mu) {
    trgswEncZeroNtt(trgsw, trgswDft, param, trgswKey);
    trgswAddIntegerNtt(trgswDft, trgsw, mu, param);
}

// trgsw(0)
void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            Trlwe& trlweSample = trgsw.trlweSamples[lvl][row];
            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[lvl][row];
            initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, param.lweStdDev); // init b = 0 + e
            applyNtt(trlweDftSample.b, trlweSample.b);
            for (auto col = 0; col < param.k; col++) {
                initCoeffsViaUniformDistribution(trlweSample.a[col].coeffs); // init a
                applyNtt(trlweDftSample.a[col], trlweSample.a[col]);
                applyNtt(trgswKey.trlweKey.sDft[col], trgswKey.trlweKey.s[col]);
                modularAccumulate(trlweDftSample.b.coeffs, trlweDftSample.a[col].coeffs, trgswKey.trlweKey.sDft[col].coeffs);
            }
            applyIntt(trlweSample.b, trlweDftSample.b);
        }
    }
}

// output += mu * G^T
void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, Integer mu, const YatfheParameters& param) {
    // add the diagonal matrix (mu * G^T)_ijk to the output
    //       ( 1/B^l                         )
    //      .                              . .
    //    .                              .   .
    //  ( 1/B^2                        )     .
    // ( 1/B                         )       .
    // (     1/B                     )       .
    // (          .                  )       .
    // (              .              )     .
    // (                  .          )   .
    // (                      .      ) .
    // (                         1/B )
    // ( a_0  a_1          a_k-1  b  )

    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            auto decomposedMu = mu << (param.torusBits - param.radixBits * (lvl + 1));

            // add to a_lii
            if (row < param.k) {
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
//            LagrangePolynomial tmp(param.N);
//            generateLagrangePolynomialWithValueAt(tmp, mu * g[layer], 0);
//            lagrangePolynomialAccumulate(output.trlweDftSamples[layer][row].b, tmp);
        }
    }
}






