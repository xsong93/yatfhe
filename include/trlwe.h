//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include <vector>
#include "torus.h"
#include "polynomial.h"

struct Trlwe {
    TorusPolynomial* a{ new TorusPolynomial }; // k
    TorusPolynomial* b{ new TorusPolynomial }; // 1
//    int k;
};

struct TrlweKey {
//    IntPolynomial* s{ new IntPolynomial[k] }; // k
    std::vector<IntPolynomial> s; // k
    int k{};
//    double sigma;
};

struct TrlweDft{
    LagrangePolynomial* a{ new LagrangePolynomial };
    LagrangePolynomial* b{ new LagrangePolynomial };
//    int k;
};

void trlweNewBinaryKey(TrlweKey& trlweKey, int N, int k);

void trlweInitKey(TrlweKey& key, int N, int k);

void initTrlweSample(Trlwe& trlwe, int k, int N);

void initTrlweDftSample(TrlweDft& trlweDft, int k, int N);

void trlweKeyGen(TrlweKey& key, int N, int k);

void deleteRlweKey(TrlweKey& key);

void deleteRlweSample(Trlwe& sample);

void deleteRlweDftSample(TrlweDft& sample);

#endif //HLS_YATFHE_TRLWE_H
