//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include <vector>
#include "tlwe.h"
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

struct TrlweDft{
    std::vector<LagrangePolynomial> a {}; // k + 1
    LagrangePolynomial b; // 1
    int k;

    TrlweDft(int k, int N) :
            a(k + 1, LagrangePolynomial(N)),
            b(LagrangePolynomial(N)),
            k(k) {};
};

struct TrlweKey {
    std::vector<IntPolynomial> s {}; // k
    std::vector<LagrangePolynomial> sDft {}; // k
    int k {};
//    double sigma;

    TrlweKey(int k, int N):
        k(k),
        s(k, TorusPolynomial(N)),
        sDft(k, LagrangePolynomial(N)) {};
};

//void newBinaryTrlweKey(TrlweKey& trlweKey, const YatfheParameters& param);
//
//void initTrlweKey(TrlweKey& key, int N, int k);
//
//void initTrlweSample(Trlwe& trlwe, int k, int N);
//
//void initTrlweDftSample(TrlweDft& trlweDft, int k, int N);

void trlweKeyGen(TrlweKey& key, int N, int k);

void trlweAccumulate(Trlwe& res, const Trlwe& accum);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

//void deleteRlweKey(TrlweKey& key);
//
//void deleteRlweSample(Trlwe& sample);
//
//void deleteRlweDftSample(TrlweDft& sample);

#endif //HLS_YATFHE_TRLWE_H
