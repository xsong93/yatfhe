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
    cout << "TrlweKey: ";
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < N; j++) {
            key.s[i].coeffs[j] = distribution(rng);
            cout << j << ":" <<key.s[i].coeffs[j] <<" ";
        }
    }
    cout << endl;
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

void genNoiselessTrlweSample(Trlwe& acc, const Torus msg, const NegaCyclicTlwe& negaCyclicInput, const YatfheParameters& param) {
    const auto& barb = negaCyclicInput.b;
    if (barb == 0) {
        std::fill(acc.b.coeffs.begin(), acc.b.coeffs.end(), msg);
        return;
    }
    const int N = param.N;
    const int N2 = negaCyclicInput.N2;
    const int rot = N2 - barb;
    const int rotTrue = (rot < N) ? rot : rot - N;
    for (int i = 0; i < N; i++) {
        acc.b.coeffs[i] = (i < rotTrue) ? -msg : msg;
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