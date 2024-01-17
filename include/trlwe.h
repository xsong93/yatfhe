//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include <vector>
#include "yatfhe_parameters.h"
#include "torus.h"
#include "polynomial.h"

struct Trlwe {
    std::vector<TorusPolynomial> a {}; // k
    TorusPolynomial b {}; // 1
    int k;

    Trlwe(int k, int N) :
        a(k, TorusPolynomial(N)),
        b(TorusPolynomial(N)),
        k(k) {};
};

struct TrlweKey {
    std::vector<IntPolynomial> s {}; // k
    int k {};
//    double sigma;

    TrlweKey() : k(0), s(0) {};

    explicit TrlweKey(int k, int N):
        k(k),
        s(k, TorusPolynomial(N)) {};
};

struct TrlweDft{
    std::vector<LagrangePolynomial> a {}; // k
    LagrangePolynomial b; // 1
    int k;

    TrlweDft(int k, int N) :
        a(k, LagrangePolynomial(N)),
        b(LagrangePolynomial(N)),
        k(k) {};
};

void newBinaryTrlweKey(TrlweKey& trlweKey, const YatfheParameters& yatfheParameters);

void initTrlweKey(TrlweKey& key, const int N, const int k);

void initTrlweSample(Trlwe& trlwe, int k, int N);

void initTrlweDftSample(TrlweDft& trlweDft, int k, int N);

void trlweKeyGen(TrlweKey& key, int N, int k);

void deleteRlweKey(TrlweKey& key);

void deleteRlweSample(Trlwe& sample);

void deleteRlweDftSample(TrlweDft& sample);

#endif //HLS_YATFHE_TRLWE_H
