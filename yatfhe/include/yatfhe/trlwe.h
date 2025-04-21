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
    int N;
    int k;

    Trlwe(int k, int N) :
        a(k, TorusPolynomial(N)),
        b(TorusPolynomial(N)),
        N(N),
        k(k) {};

    Trlwe(int k, int N, Torus val) :
            a(k, TorusPolynomial(N, val)),
            b(TorusPolynomial(N, val)),
            N(N),
            k(k) {};

    Trlwe(int k, int N, Torus valA, Torus valB) :
            a(k, TorusPolynomial(N, valA)),
            b(TorusPolynomial(N, valB)),
            N(N),
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

struct Trlwe8D {
    std::vector<Int8PolynomialD> a; // k
    Int8PolynomialD b; // 1
    int k;

    Trlwe8D(int k, int N, int d) :
            a(k, Int8PolynomialD(N, d)),
            b(Int8PolynomialD(N, d)),
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

struct TrlweDft24{
    std::vector<Ntt24Polynomial> a; // k
    Ntt24Polynomial b; // 1
    int k;

    TrlweDft24(int k, int N) :
            a(k, Ntt24Polynomial(N)),
            b(Ntt24Polynomial(N)),
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
            s(param.k, TorusPolynomial(param.N)),
            sDft(param.k, LagrangePolynomial(param.N)),
            k(param.k),
            N(param.N),
            sigma(param.rlweStdDev) {};

    TrlweKey(int k, int N, double sigma):
        s(k, TorusPolynomial(N)),
        sDft(k, LagrangePolynomial(N)),
        k(k),
        N(N),
        sigma(sigma) {};
};

template<typename TrlweType>
void trlweAdd(TrlweType& output, const TrlweType& input1, const TrlweType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        polynomialAddT32(output.a[i], input1.a[i], input2.a[i]);
    }
    polynomialAddT32(output.b, input1.b, input2.b);
}

template<typename TrlweType>
void trlweSub(TrlweType& output, const TrlweType& input1, const TrlweType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        polynomialSubT32(output.a[i], input1.a[i], input2.a[i]);
    }
    polynomialSubT32(output.b, input1.b, input2.b);
}

template<typename TrlweDftType>
void trlweAddNtt(TrlweDftType& output, const TrlweDftType& input1, const TrlweDftType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        lagrangePolynomialAdd(output.a[i], input1.a[i], input2.a[i]);
    }
    lagrangePolynomialAdd(output.b, input1.b, input2.b);
}

template<typename TrlweDftType>
void trlweSubNtt(TrlweDftType& output, const TrlweDftType& input1, const TrlweDftType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        lagrangePolynomialSub(output.a[i], input1.a[i], input2.a[i]);
    }
    lagrangePolynomialSub(output.b, input1.b, input2.b);
}

template <typename T>
void trlweSetZero(std::vector<T>& a, T& b) {
    std::fill(a.begin(), a.end(), T(b.N, 0));
    std::fill(b.coeffs.begin(), b.coeffs.end(), 0);
}

/**
 * accum.a += tlwe.a, accum.b += tlwe.b
 * */
template<typename TrlweType>
void trlweAccumulateI32(TrlweType& accum, const TrlweType& tlwe) {
    for (auto i = 0; i < accum.a.size(); i++) {
        polynomialAccumulateI32(accum.a[i], tlwe.a[i]);
    }
    polynomialAccumulateI32(accum.b, tlwe.b);
}

template<typename TrlweType>
void trlweAccumulateT32(TrlweType& accum, const TrlweType& tlwe) {
    for (auto i = 0; i < accum.a.size(); i++) {
        polynomialAccumulateT32(accum.a[i], tlwe.a[i]);
    }
    polynomialAccumulateT32(accum.b, tlwe.b);
}

template<typename TrlweType, typename U>
void trlweAccumulateModP(TrlweType& accum, const TrlweType& tlwe, const U p) {
    for (auto i = 0; i < accum.a.size(); i++) {
        polynomialAccumulateModP(accum.a[i], tlwe.a[i], p);
    }
    polynomialAccumulateModP(accum.b, tlwe.b, p);
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweCRTDecomp(std::vector<TrlweTypeA>& out, const TrlweTypeB& in, const YatfheParameters& param) {
    for (size_t d = 0; d < param.d; d++) {
        auto qd = param.qd[d];
        auto& outA = out[d].a;
        auto& outB = out[d].b;
        auto& inA = in.a;
        auto& inB = in.b;
        for (size_t k = 0; k < param.k; k++) {
            auto& coeffOutA = outA[k].coeffs;
            auto& coeffInA = inA[k].coeffs;
            for (size_t j = 0; j < param.N; j++) {
                coeffOutA[j] = longModP(coeffInA[j], qd);
            }
        }
        auto& coeffOutB = outB.coeffs;
        auto& coeffInB = inB.coeffs;
        for (size_t j = 0; j < param.N; j++) {
            coeffOutB[j] = longModP(coeffInB[j], qd);
        }
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweCRTDecompNO(TrlweTypeA& out, const TrlweTypeB& in, const YatfheParameters& param) {
    auto& outA = out.a;
    auto& outB = out.b;
    auto& inA = in.a;
    auto& inB = in.b;
    for (size_t k = 0; k < param.k; k++) {
        auto& coeffOutA = outA[k].coeffs;
        auto& coeffInA = inA[k].coeffs;
        for (size_t j = 0; j < param.N; j++) {
            auto& coeffOutDA = coeffOutA[j];
            auto valA = coeffInA[j];
            for (size_t d = 0; d < param.d; d++) {
                coeffOutDA[d] = longModP(valA, param.qd[d]);
            }
        }
    }
    auto& coeffOutB = outB.coeffs;
    auto& coeffInB = inB.coeffs;
    for (size_t j = 0; j < param.N; j++) {
        auto& coeffOutDB = coeffOutB[j];
        auto valB = coeffInB[j];
        for (size_t d = 0; d < param.d; d++) {
            coeffOutDB[d] = longModP(valB, param.qd[d]);
        }
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweMCRTDecomp(std::vector<TrlweTypeA>& out, const TrlweTypeB& in, const YatfheParameters& param) {
    for (size_t d = 0; d < param.d; d++) {
        int64_t taoU = param.taoU[d];
        auto qd = param.qd[d];
        auto& outA = out[d].a;
        auto& outB = out[d].b;
        auto& inA = in.a;
        auto& inB = in.b;
        for (size_t k = 0; k < param.k; k++) {
            auto& coeffOutA = outA[k].coeffs;
            auto& coeffInA = inA[k].coeffs;
            for (size_t j = 0; j < param.N; j++) {
                coeffOutA[j] = static_cast<int8_t>(longModP(taoU * coeffInA[j], qd));
            }
        }
        auto& coeffOutB = outB.coeffs;
        auto& coeffInB = inB.coeffs;
        for (size_t j = 0; j < param.N; j++) {
            coeffOutB[j] = static_cast<int8_t>(longModP(taoU * coeffInB[j], qd));
        }
    }
}

template<typename TrlweType>
void trlweMCRTToCRT(std::vector<TrlweType>& trlwe, const YatfheParameters& param) {
    auto dh = param.dh;
    for (size_t d = 0; d < param.dl; d++) {
        auto ql = param.ql[d];
        int64_t taoUInv = param.taoUInv[dh + d];
        auto& accA = trlwe[dh + d].a;
        auto& accB = trlwe[dh + d].b;
        for (size_t k = 0; k < param.k; k++) {
            auto& coeffA = accA[k].coeffs;
            for (size_t j = 0; j < param.N; j++) {
                auto aCopy = coeffA[j];
                coeffA[j] = static_cast<int8_t>(longModP(taoUInv * aCopy, ql));
            }
        }
        for (size_t j = 0; j < param.N; j++) {
            auto bCopy = accB.coeffs[j];
            accB.coeffs[j] = static_cast<int8_t>(longModP(taoUInv * bCopy, ql));
        }
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweApproxCRTRecomp(TrlweTypeA& out, std::vector<TrlweTypeB>& inMCRT, const YatfheParameters& param) {
    auto& outA = out.a;
    auto qCRT = param.qCRT;
    for (size_t k = 0; k < param.k; k++) {
        auto& coeffA = outA[k].coeffs;
        for (size_t j = 0; j < param.N; j++) {
            long tmpA = 0;
            for (size_t d = 0; d < param.dh; d++) {
                tmpA += inMCRT[d].a[k].coeffs[j] * param.w[d];
            }
            coeffA[j] = static_cast<Torus>(longModP(tmpA, qCRT));
        }
    }
    auto& coeffB = out.b.coeffs;
    for (size_t j = 0; j < param.N; j++) {
        long tmpB = 0;
        for (size_t d = 0; d < param.dh; d++) {
            tmpB += inMCRT[d].b.coeffs[j] * param.w[d];
        }
        coeffB[j] = static_cast<Torus>(longModP(tmpB, qCRT));
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweCRTRecomp(TrlweTypeA& out, std::vector<TrlweTypeB>& inCRT, const YatfheParameters& param) {
    auto& outA = out.a;
    auto qCRT = param.qCRT;
    for (size_t k = 0; k < param.k; k++) {
        auto& coeffA = outA[k].coeffs;
        for (size_t j = 0; j < param.N; j++) {
            long tmpA = 0;
            for (size_t d = 0; d < param.d; d++) {
                tmpA += inCRT[d].a[k].coeffs[j] * param.z[d];
            }
            coeffA[j] = static_cast<Torus>(longModP(tmpA, qCRT));
        }
    }
    auto& coeffB = out.b.coeffs;
    for (size_t j = 0; j < param.N; j++) {
        long tmpB = 0;
        for (size_t d = 0; d < param.d; d++) {
            tmpB += inCRT[d].b.coeffs[j] * param.z[d];
        }
        coeffB[j] = static_cast<Torus>(longModP(tmpB, qCRT));
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void trlweCRTRecompNO(TrlweTypeA& out, TrlweTypeB& inCRT, const YatfheParameters& param) {
    auto& outA = out.a;
    auto& inA = inCRT.a;
    auto qCRT = param.qCRT;
    for (size_t k = 0; k < param.k; k++) {
        auto& coeffOutA = outA[k].coeffs;
        auto& coeffInA = inA[k].coeffs;
        for (size_t j = 0; j < param.N; j++) {
            long tmpA = 0;
            auto& coeffInDA = coeffInA[j];
            for (size_t d = 0; d < param.d; d++) {
                tmpA += coeffInDA[d] * param.z[d];
            }
            coeffOutA[j] = static_cast<Torus>(longModP(tmpA, qCRT));
        }
    }
    auto& coeffOutB = out.b.coeffs;
    auto& coeffInB = inCRT.b.coeffs;
    for (size_t j = 0; j < param.N; j++) {
        long tmpB = 0;
        auto& coeffInDB = coeffInB[j];
        for (size_t d = 0; d < param.d; d++) {
            tmpB += coeffInDB[d] * param.z[d];
        }
        coeffOutB[j] = static_cast<Torus>(longModP(tmpB, qCRT));
    }
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

//void trlweAccumulateI32(Trlwe& accum, const Trlwe& tlwe);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param);

void gadgetDecomposeTrlweNtt(DecomposedTrlweDft& output, const TrlweDft& input, const YatfheParameters& param);

void recomposeTrlwe(Trlwe& output, const DecomposedTrlwe& input, const YatfheParameters& param);

void recomposeTrlweNtt(TrlweDft& output, const DecomposedTrlweDft& input, const YatfheParameters& param);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey);

void trlweRotate(Trlwe& res, const Trlwe& input, int a);

void trlweRotateNtt(TrlweDft& res, const TrlweDft& input, const int r);

void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, int a);

void trlweRotateMinusOne8(Trlwe8& res, const Trlwe8& input, int a, int modP);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

#endif //HLS_YATFHE_TRLWE_H
