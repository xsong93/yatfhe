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
#include "yatfhe/gadget_decomposition.h"

struct Rlwe {
    std::vector<IntPolynomial> a {}; // k
    IntPolynomial b; // 1
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

    explicit DecomposedTrlwe(const YatfheParameters& param) :
            l(param.l),
            lDft(param.l * (param.dftBits / param.torusBits)),
            rlwes(param.l,  Rlwe(param.k, param.N)),
            rlweDfts(param.l * (param.dftBits / param.torusBits), TrlweDft(param.k, param.N)) {};
};

struct TrlweKey {
    std::vector<BinPolynomial> s; // k
    std::vector<LagrangePolynomial> sDft; // k
    int k;
    int N;
    double sigma {};

    explicit TrlweKey(const YatfheParameters& param):
            k(param.k),
            N(param.N),
            sigma(param.rlweStdDev),
            s(param.k, TorusPolynomial(param.N)),
            sDft(param.k, LagrangePolynomial(param.N)) {};

    TrlweKey(int k, int N, double sigma):
        k(k),
        N(N),
        sigma(sigma),
        s(k, TorusPolynomial(N)),
        sDft(k, LagrangePolynomial(N)) {};
};

template <typename T>
void trlweSetZero(std::vector<T>& a, T& b) {
    std::fill(a.begin(), a.end(), T(b.N, 0));
    std::fill(b.coeffs.begin(), b.coeffs.end(), 0);
}

void trlweKeyGen(TrlweKey& key);

void symEncTrlweSingleSample(Trlwe& trlwe, const TrlweKey& key, Torus mu);

void symEncTrlweMultiSample(Trlwe& trlwe, const TrlweKey& key, const std::vector<Torus>& mu);

void symEncTrlweSingleSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, Torus mu);

void symEncTrlweMultiSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const std::vector<Torus>& mu);

void symDecTrlweToDouble(DoublePolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweToTorus(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweToInt(IntPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweWoRounding(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key);

void symDecTrlweNtt(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRoundingNtt(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key);

void trlweAdd(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweSub(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweAddNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweSubNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param);

void gadgetDecomposeTrlweNtt(DecomposedTrlwe& output, const TrlweDft& input, const YatfheParameters& param);

void recomposeTrlwe(Trlwe& output, const DecomposedTrlwe& input, const YatfheParameters& param);

void recomposeTrlweNtt(TrlweDft& output, const DecomposedTrlwe& input, const YatfheParameters& param);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

#endif //HLS_YATFHE_TRLWE_H
