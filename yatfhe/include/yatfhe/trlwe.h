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
    std::vector<TrlweDft> rlweDfts; // 2l
    int l;
    int lDft;

    explicit DecomposedTrlwe(YatfheParameters param) :
            l(param.l),
            lDft(param.l * (param.dftBits / param.torusBits)),
            rlwes(param.l,  Rlwe(param.k, param.N)),
            rlweDfts(param.l * (param.dftBits / param.torusBits), TrlweDft(param.k, param.N)) {};
};

struct Tglev {
    std::vector<Trlwe> trlwes; // l
    std::vector<TrlweDft> trlweDfts; // 2l
    int l;
    int lDft;

    explicit Tglev(YatfheParameters param) :
            l(param.l),
            lDft(param.l * (param.dftBits / param.torusBits)),
            trlwes(param.l,  Trlwe(param.k, param.N)),
            trlweDfts(param.l * (param.dftBits / param.torusBits), TrlweDft(param.k, param.N)) {};
};

struct TrlweKey {
    std::vector<BinPolynomial> s; // k
    std::vector<LagrangePolynomial> sDft; // k
    int k;
    int N;

    TrlweKey(int k, int N):
        k(k),
        N(N),
        s(k, TorusPolynomial(N)),
        sDft(k, LagrangePolynomial(N)) {};
};

template <typename T>
void trlweSetZero(std::vector<T>& a, T& b) {
    std::fill(a.begin(), a.end(), T(b.N, 0));
    std::fill(b.coeffs.begin(), b.coeffs.end(), 0);
}

void trlweKeyGen(TrlweKey& key);

void symEncTrlweSingleSample(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, Torus mu, double sigma);

void symEncTrlweMultiSample(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const std::vector<Torus>& mu, double sigma);

void symDecTrlwe(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRounding(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key);

void trlweAdd(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweSub(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweAddNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweSubNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

#endif //HLS_YATFHE_TRLWE_H
