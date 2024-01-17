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
//    TorusPolynomial* a{ new TorusPolynomial }; // k
    std::vector<TorusPolynomial> a {}; // k
    TorusPolynomial b {}; // 1
//    int k;
};

struct TrlweKey {
//    IntPolynomial* bskDft{ new IntPolynomial[k] }; // k
    std::vector<IntPolynomial> s {}; // k
    int k {};
//    double sigma;
};

struct TrlweDft{
//    LagrangePolynomial* a{ new LagrangePolynomial };
//    LagrangePolynomial* b{ new LagrangePolynomial };
    std::vector<LagrangePolynomial> a {}; // k
    LagrangePolynomial b {}; // 1
//    int k;
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
