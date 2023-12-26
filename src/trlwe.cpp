//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "trlwe.h"
#include "polynomial.h"
#include "numeric_functions.h"

using namespace std;

TrlweKey *trlweNewBinaryKey(const int N, const int k) {
    TrlweKey* trlweKey{new TrlweKey};
    trlweKey->k = k;
//    TrlweKey* trlweKey = trlweInitKey(N, k, sigma);
    trlweInitKey(trlweKey, N, k);
    trlweKeyGen(trlweKey, N, k);
    return trlweKey;
}

void trlweInitKey(TrlweKey* key, const int N, const int k) {
//    key->sigma = sigma;
    key->s = new IntPolynomial;
    for (int i = 0; i < k; i++) {
//        key->s[i] = *initTorusPolynomial(N);
        initTorusPolynomial(&key->s[i], N);
//        std::cout << key->s[i].N << " ";
    }
//    std::cout <<endl;
}

//TrlweKey* trlweInitKey(int N, int k, double sigma) {
//    TrlweKey* key{new TrlweKey};
//    key->k = k;
//    key->sigma = sigma;
//    key->s = new IntPolynomial;
//    for (int i = 0; i < k; i++) {
////        key->s[i] = *initTorusPolynomial(N);
//        initTorusPolynomial(&key->s[i], N);
////        std::cout << key->s[i].N << " ";
//    }
////    std::cout <<endl;
//    return key;
//}

void trlweKeyGen(TrlweKey* key, const int N, const int k) {
    uniform_int_distribution<int> distribution(0, 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < N; j++) {
            key->s[i].coeffs[j] = distribution(generator);
//            cout << key->s[i].coeffs[j] <<" ";
        }
    }
//    cout <<endl;
}

//Trlwe* initTrlweSample(int k, int N) {
//    Trlwe* res{new Trlwe};
//    res->a = new TorusPolynomial;
//    for (int i = 0; i < k; i++) {
//        res->a[i] = *initTorusPolynomial(N);
//    }
//    res->b = initTorusPolynomial(N);
//    res->k = k;
//    return res;
//}

void initTrlweSample(Trlwe* trlwe, const int k, const int N) {
    trlwe->a = new TorusPolynomial;
    for (int i = 0; i < k; i++) {
//        trlwe->a[i] = *initTorusPolynomial(N);
        initTorusPolynomial(&trlwe->a[i], N);
    }
//    trlwe->b = initTorusPolynomial(N);
    initTorusPolynomial(trlwe->b, N);
//    trlwe->k = k;
}

//TrlweDft* initTrlweDftSample(int k, int N) {
//    TrlweDft* res{new TrlweDft};
//    res->a = new LagrangePolynomial;
//    for (int i = 0; i < k; i++) {
//        res->a[i] = *initLagrangePolynomial(N);
//    }
//    res->b = initLagrangePolynomial(N);
//    res->k = k;
//    return res;
//}

void initTrlweDftSample(TrlweDft* trlweDft, const int k, const int N) {
    trlweDft->a = new LagrangePolynomial;
    for (int i = 0; i < k; i++) {
//        trlweDft->a[i] = *initLagrangePolynomial(N);
        initLagrangePolynomial(&trlweDft->a[i], N);
    }
//    trlweDft->b = initLagrangePolynomial(N);
    initLagrangePolynomial(trlweDft->b, N);
//    trlweDft->k = k;
}

void deleteRlweKey(TrlweKey* key) {
    delete key->s->coeffs;
    key->s->coeffs = nullptr;
    delete key->s;
    key->s = nullptr;
    delete key;
    key = nullptr;
}