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

    Trlwe(int k, int N, int val) :
            a(k, TorusPolynomial(N, val)),
            b(TorusPolynomial(N, val)),
            k(k) {};
};

struct Trlwe8 {
    std::vector<Int8Polynomial> a; // k
    Int8Polynomial b; // 1
    int k;

    Trlwe8(int k, int N) :
            a(k, Int8Polynomial(N)),
            b(Int8Polynomial(N)),
            k(k) {};

    Trlwe8(int k, int N, int8_t val) :
            a(k, Int8Polynomial(N, val)),
            b(Int8Polynomial(N, val)),
            k(k) {};
};

struct Trlwe16 {
    std::vector<Int16Polynomial> a; // k
    Int16Polynomial b; // 1
    int k;

    Trlwe16(int k, int N) :
            a(k, Int16Polynomial(N)),
            b(Int16Polynomial(N)),
            k(k) {};

    Trlwe16(int k, int N, int16_t val) :
            a(k, Int16Polynomial(N, val)),
            b(Int16Polynomial(N, val)),
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

struct TrlweDft14{
    std::vector<Ntt14Polynomial> a; // k
    Ntt14Polynomial b; // 1
    int k;

    TrlweDft14(int k, int N) :
            a(k, Ntt14Polynomial(N)),
            b(Ntt14Polynomial(N)),
            k(k) {};
};

struct TrlweDft16{
    std::vector<Ntt16Polynomial> a; // k
    Ntt16Polynomial b; // 1
    int k;

    TrlweDft16(int k, int N) :
            a(k, Ntt16Polynomial(N)),
            b(Ntt16Polynomial(N)),
            k(k) {};
};

struct DecomposedTrlwe {
    std::vector<Rlwe> rlwes; // l
    int l;

    explicit DecomposedTrlwe(const YatfheParameters& param) :
            l(param.l),
            rlwes(param.l,  Rlwe(param.k, param.N)) {};
    DecomposedTrlwe(const YatfheParameters& param, int l) :
            l(l),
            rlwes(l,  Rlwe(param.k, param.N)) {};
};

struct DecomposedTrlweDft {
    std::vector<TrlweDft> rlweDfts; // 2l, l
    int l;

    DecomposedTrlweDft(const YatfheParameters& param, int l) :
            l(l),
            rlweDfts(l, TrlweDft(param.k, param.N)) {};
};

struct DecomposedTrlweDft14 {
    std::vector<TrlweDft14> rlweDfts; // 2l, l
    int l;

    DecomposedTrlweDft14(const YatfheParameters& param, int l) :
            l(l),
            rlweDfts(l, TrlweDft14(param.k, param.N)) {};
};

struct DecomposedTrlweDft16 {
    std::vector<TrlweDft16> rlweDfts; // 2l, l
    int l;

    DecomposedTrlweDft16(const YatfheParameters& param, int l) :
            l(l),
            rlweDfts(l, TrlweDft16(param.k, param.N)) {};
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

void symDecTrlweToIntNtt(IntPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRounding(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key);

void symDecTrlweNtt(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRoundingNtt(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key);

void trlweAdd(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweSub(Trlwe& output, const Trlwe& input1, const Trlwe& input2);

void trlweAddNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweSubNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2);

void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param);

void gadgetDecomposeTrlweNtt(DecomposedTrlweDft& output, const TrlweDft& input, const YatfheParameters& param);

void recomposeTrlwe(Trlwe& output, const DecomposedTrlwe& input, const YatfheParameters& param);

void recomposeTrlweNtt(TrlweDft& output, const DecomposedTrlweDft& input, const YatfheParameters& param);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey);

void trlweRotate(Trlwe& res, const Trlwe& input, int a);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void trlweRotateMinusOne8(Trlwe8& res, const Trlwe8& input, int a);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

#endif //HLS_YATFHE_TRLWE_H
