//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include "torus.h"
#include "polynomial.h"

struct Trlwe {
    TorusPolynomial* a; // k
    TorusPolynomial* b; // 1
//    int k;
};

struct TrlweKey {
    IntPolynomial* s;
    int k;
//    double sigma;
};

struct TrlweDft{
    LagrangePolynomial* a;
    LagrangePolynomial* b;
//    int k;
};

TrlweKey *trlweNewBinaryKey(int N, int k);

void trlweInitKey(TrlweKey *key, int N, int k);

//TrlweKey* trlweInitKey(int N, int k, double sigma);

//Trlwe* initTrlweSample(int k, int N);

void initTrlweSample(Trlwe* trlwe, int k, int N);

TrlweDft* initTrlweDftSample(int k, int N);

void initTrlweDftSample(TrlweDft* trlweDft, int k, int N);

void trlweKeyGen(TrlweKey* key, int N, int k);

void deleteRlweKey(TrlweKey* key);

#endif //HLS_YATFHE_TRLWE_H
