//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include <vector>
#include "yatfhe/tlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"

struct Rlwe {
    std::vector<IntPolynomial> a {}; // k
    IntPolynomial b {}; // 1
    int k;

    Rlwe(int k, int N) :
            a(k, IntPolynomial(N)),
            b(IntPolynomial(N)),
            k(k) {};
};

struct Trlwe {
    std::vector<TorusPolynomial> a; // k
    TorusPolynomial b; // 1
    int k;

    Trlwe(int k, int N) :
        a(k, TorusPolynomial(N)),
        b(TorusPolynomial(N)),
        k(k) {};
};

struct TrlweDft{
    std::vector<LagrangePolynomial> a; // k
    LagrangePolynomial b; // 1
    int k;

    TrlweDft(int k, int N) :
            a(k, LagrangePolynomial(N)),
            b(LagrangePolynomial(N)),
            k(k) {};
};


struct DecomposedTrlwe {
    std::vector<Rlwe> rlwes; // l
    std::vector<TrlweDft> rlweDfts; // l
    int l;

    DecomposedTrlwe(int l, int k, int N) :
            l(l),
            rlwes(l,  Rlwe(k, N)),
            rlweDfts(l, TrlweDft(k, N)) {};
};

struct TrlweKey {
    std::vector<IntPolynomial> s; // k
    std::vector<LagrangePolynomial> sDft; // k
    int k;
//    double sigma;

    TrlweKey(int k, int N):
        k(k),
        s(k, TorusPolynomial(N)),
        sDft(k, LagrangePolynomial(N)) {};
};

void trlweKeyGen(TrlweKey& key, int N, int k);

void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

void gadgetDecomposition(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_TRLWE_H
