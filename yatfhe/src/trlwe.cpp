//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "yatfhe_parameters.h"
#include "tlwe.h"
#include "trlwe.h"
#include "polynomial.h"
#include "numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

//void newBinaryTrlweKey(TrlweKey& trlweKey, const YatfheParameters& param) {
//    initTrlweKey(trlweKey, param.N, param.k);
//    trlweKeyGen(trlweKey, param.N, param.k);
//}

//void initTrlweKey(TrlweKey& key, const int N, const int k) {
//    for (int i = 0; i < k; i++) {
//        TorusPolynomial torusPolynomial(N);
//        key.s.push_back(torusPolynomial);
//    }
//}

void trlweKeyGen(TrlweKey& key, const int N, const int k) {
    uniform_int_distribution<int> distribution(0, 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < N; j++) {
            key.s[i].coeffs[j] = distribution(rng);
        }
    }
//    printPolyVec(key.s, "TrlweKey");
}

//void initTrlweSample(Trlwe& trlwe, const int k, const int N) {
////    trlwe.a.resize(k);
//    for (int i = 0; i < k; i++) {
//        initTorusPolynomial(trlwe.a[i], N);
//    }
//    initTorusPolynomial(trlwe.b, N);
//}
//
//void initTrlweDftSample(TrlweDft& trlweDft, const int k, const int N) {
////    trlweDft.a.resize(k);
//    for (int i = 0; i < k; i++) {
//        initLagrangePolynomial(trlweDft.a[i], N);
//    }
//    initLagrangePolynomial(trlweDft.b, N);
//}

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
 * res.a += accum.a
 * */
void trlweAccumulate(Trlwe& res, const Trlwe& accum) {
    const auto k = res.k;
    for (auto i = 0; i < k + 1; i++) {
        polynomialAccumulate(res.a[i], accum.a[i]);
    }
//    ploynomialAccumulate(res.b, accum.b);
}

// out = (a', b[index])
void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, const int index) {
    const auto N = in.b.N;
    const auto size = in.a.size();
    for (auto i = 0; i < size; i++) {
        for (auto j = 0; j < N; j++) {
            out.a[i * N + j] = (j <= index) ? (in.a[i].coeffs[index - j]) : (-in.a[i].coeffs[N + index - j]);
        }
    }
    out.b = in.b.coeffs[index];
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

//void deleteRlweKey(TrlweKey& key) {
////    const int k = key.k;
////    for (int i = 0; i < k; i++) {
////        delete[] key.bskDft[i].coeffs;
////        key.bskDft[i].coeffs = nullptr;
////    }
////    delete[] key.bskDft;
////    key.bskDft = nullptr;
//}
//
//void deleteRlweSample(Trlwe& sample) {
////    deletePolynomial(*sample.a);
////    deletePolynomial(*sample.b);
//}
//
//void deleteRlweDftSample(TrlweDft& sample) {
////    deletePolynomial(*sample.a);
////    deletePolynomial(*sample.b);
//}