//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yautil/tool.h"

using namespace std;

void initTrlweSingleSample(Trlwe& trlwe, const Torus mu, double sigma) {
    initCoeffsWithGaussianNoiseSingleSample(trlwe.b.coeffs, mu, sigma);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs);
    }
}

void initTrlweMultiSample(Trlwe& trlwe, const vector<Torus>& mu, double sigma) {
    initCoeffsWithGaussianNoiseMultiSample(trlwe.b.coeffs, mu, sigma);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs);
    }
}

void symEncTrlwe(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key) {
    applyNttForAB(trlweDft, trlwe);
    calModularInnerProductNtt(trlweDft.b, trlweDft.a, key.sDft);
    applyIntt(trlwe.b, trlweDft.b);
}

void trlweKeyGen(TrlweKey& key) {
    for (int i = 0; i < key.k; i++) {
        for (int j = 0; j < key.N; j++) {
            key.s[i].coeffs[j] = binaryDistrib(rng);
        }
        applyNtt(key.sDft[i], key.s[i]);
    }
//    printPolyVec(key.s, "TrlweKey");
}

void symEncTrlweSingleSample(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const Torus mu, double sigma) {
    initTrlweSingleSample(trlwe, mu, sigma);
    symEncTrlwe(trlwe, trlweDft, key);
}

void symEncTrlweMultiSample(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const vector<Torus>& mu, double sigma) {
    initTrlweMultiSample(trlwe, mu, sigma);
    symEncTrlwe(trlwe, trlweDft, key);
}

void symDecTrlwe(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    LagrangePolynomial innerProduct {trlweDft.b.N};
    LagrangePolynomial res {trlweDft.b.N};
    calModularInnerProductNtt(innerProduct, trlweDft.a, key.sDft);
    lagrangePolynomialSub(res, trlweDft.b, innerProduct);
    applyIntt(tmp, res);
    torusPolyToDoublePoly(output, tmp);
    roundErrorPoly(output, torusBase);
}

void symDecTrlweWoRounding(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    DoublePolynomial tmp {output.N};
    LagrangePolynomial innerProduct {trlweDft.b.N};
    LagrangePolynomial res {trlweDft.b.N};
    calModularInnerProductNtt(innerProduct, trlweDft.a, key.sDft);
    lagrangePolynomialSub(res, trlweDft.b, innerProduct);
    applyIntt(output, res);
}

// Trlwe: (X^-b) * (0,...,0,v)
void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput) {
    const auto barb = scaledInput.b;
    torusPolynomialRotate(accum.b, -barb, v);

    std::vector<double> t(accum.b.coeffs.size());
    for (int i = 0; i < accum.b.N; i++) {
        t[i] = torus32ToDouble(accum.b.coeffs[i]);
    }
    printArray(t, "b:");
}

/**
 * accum.a += tlwe.a, accum.b += tlwe.b
 * */
void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe) {
    for (auto i = 0; i < accum.a.size(); i++) {
        polynomialAccumulate(accum.a[i], tlwe.a[i]);
    }
    polynomialAccumulate(accum.b, tlwe.b);
}

// out = (a', b[index])
void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, const int index) {
    const auto N = in.b.N;
    for (auto i = 0; i < in.k; i++) {
        for (auto j = 0; j < N; j++) {
            out.a[i * N + j] = (j <= index) ? (in.a[i].coeffs[index - j]) : (-in.a[i].coeffs[N + index - j]);
        }
    }
    out.b = in.b.coeffs[index];
}

// flatten a trlwe key as a tlwe key
void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey) {
    for (auto i = 0; i < trlweKey.k; i++) {
        const auto N = trlweKey.s[i].N;
        for (auto j = 0; j < N; j++) {
            tlweKey.s[i * N + j] = trlweKey.s[i].coeffs[j];
        }
    }
}

// res = X^a * input - input
void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, const int a) {
    const auto size = input.a.size();
    for (auto i = 0; i < size; i++) {
        torusPolynomialRotateMinusOne(res.a[i], a, input.a[i]);
    }
    torusPolynomialRotateMinusOne(res.b, a, input.b);
}

void copyTrlwe(Trlwe& target, const Trlwe& source, const bool copyA, const bool copyB) {
    const auto size = source.a.size();
    const auto N = source.b.N;
    for (auto j = 0; j < N; j++) {
        if (copyA) {
            for (auto i = 0; i < size; i++) {
                target.a[i].coeffs[j] = source.a[i].coeffs[j];
            }
        }
        if (copyB) {
            target.b.coeffs[j] = source.b.coeffs[j];
        }
    }
}
