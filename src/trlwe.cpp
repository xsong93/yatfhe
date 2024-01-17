//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "yatfhe_parameters.h"
#include "trlwe.h"
#include "polynomial.h"
#include "numeric_functions.h"

using namespace std;

void newBinaryTrlweKey(TrlweKey& trlweKey, const YatfheParameters& yatfheParameters) {
    trlweKey.k = yatfheParameters.k;
    initTrlweKey(trlweKey, yatfheParameters.N, yatfheParameters.k);
    trlweKeyGen(trlweKey, yatfheParameters.N, yatfheParameters.k);
}

void initTrlweKey(TrlweKey& key, const int N, const int k) {
    key.s.resize(k);
    for (int i = 0; i < k; i++) {
        initTorusPolynomial(key.s[i], N);
    }
}

void trlweKeyGen(TrlweKey& key, const int N, const int k) {
    uniform_int_distribution<int> distribution(0, 1);
    cout << "TrlweKey: ";
    for (int i = 0; i < k; i++) {
        key.s[i].coeffs.resize(N);
        for (int j = 0; j < N; j++) {
            key.s[i].coeffs[j] = distribution(rng);
            cout << j << ":" <<key.s[i].coeffs[j] <<" ";
        }
    }
    cout << endl;
}

void initTrlweSample(Trlwe& trlwe, const int k, const int N) {
    trlwe.a.resize(k);
    for (int i = 0; i < k; i++) {
        initTorusPolynomial(trlwe.a[i], N);
    }
    initTorusPolynomial(trlwe.b, N);
}

void initTrlweDftSample(TrlweDft& trlweDft, const int k, const int N) {
    trlweDft.a.resize(k);
    for (int i = 0; i < k; i++) {
        initLagrangePolynomial(trlweDft.a[i], N);
    }
    initLagrangePolynomial(trlweDft.b, N);
}

void deleteRlweKey(TrlweKey& key) {
//    const int k = key.k;
//    for (int i = 0; i < k; i++) {
//        delete[] key.bskDft[i].coeffs;
//        key.bskDft[i].coeffs = nullptr;
//    }
//    delete[] key.bskDft;
//    key.bskDft = nullptr;
}

void deleteRlweSample(Trlwe& sample) {
//    deletePolynomial(*sample.a);
//    deletePolynomial(*sample.b);
}

void deleteRlweDftSample(TrlweDft& sample) {
//    deletePolynomial(*sample.a);
//    deletePolynomial(*sample.b);
}