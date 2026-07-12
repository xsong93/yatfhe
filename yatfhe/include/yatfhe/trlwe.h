//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRLWE_H
#define HLS_YATFHE_TRLWE_H

#include <vector>

#include "ntt_hexl.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/torus.h"
#include "yatfhe/numeric.h"
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
    int N{};
    int k{};

    Trlwe() = default;

    explicit Trlwe(const YatfheParameters& p):
        a(p.k, TorusPolynomial(p.N)),
        b(TorusPolynomial(p.N)),
        N(p.N),
        k(p.k) {};

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

struct TrlweDft {
    std::vector<NttPolynomial> a;
    NttPolynomial b;
    int k{};
    bool a_initialized{};

    TrlweDft() = default;

    explicit TrlweDft(const YatfheParameters& p):
            a(p.k, NttPolynomial(p.N)),
            b(NttPolynomial(p.N)),
            k(p.k),
            a_initialized(true) {}

    TrlweDft(int k, int N) :
        a(k, NttPolynomial(N)),
        b(NttPolynomial(N)),
        k(k),
        a_initialized(true) {}

    TrlweDft(int k, int N, bool onlyB) :
        b(NttPolynomial(N)),
        k(k),
        a_initialized(!onlyB) {
        if (a_initialized) {
            a.resize(k, NttPolynomial(N));
        }
    }
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
    std::vector<Trlwe> trlwes; // l
    int l;

    explicit DecomposedTrlwe(const YatfheParameters& param) :
            l(param.l),
            trlwes(param.l,  Trlwe(param.k, param.N)) {};
    DecomposedTrlwe(const YatfheParameters& param, int l) :
            l(l),
            trlwes(l,  Trlwe(param.k, param.N)) {};
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
    std::vector<NttPolynomial> sDft; // k
    int k;
    int N;
    int errorB {};

    explicit TrlweKey(const YatfheParameters& param):
            s(param.k, BinPolynomial(param.N)),
            sDft(param.k, NttPolynomial(param.N)),
            k(param.k),
            N(param.N),
            errorB(param.rlweNoiseB) {};

    TrlweKey(int k, int N, int sigma):
        s(k, BinPolynomial(N)),
        sDft(k, NttPolynomial(N)),
        k(k),
        N(N),
        errorB(sigma) {};
};

// template<typename TrlweType>
// void addTrlwe(TrlweType& output, const TrlweType& input1, const TrlweType& input2) {
//     for (auto i = 0; i < output.a.size(); i++) {
//         addTorusPolynomial(output.a[i], input1.a[i], input2.a[i]);
//     }
//     addTorusPolynomial(output.b, input1.b, input2.b);
// }

template<typename TrlweType, typename... TrlweArgs>
void addTrlwe(TrlweType& output, const TrlweArgs&... inputs) {
    for (auto i = 0; i < output.a.size(); i++) {
        addTorusPolynomial(output.a[i], inputs.a[i]...);
    }
    addTorusPolynomial(output.b, inputs.b...);
}

template<typename TrlweType>
void subTrlwe(TrlweType& output, const TrlweType& input1, const TrlweType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        subTorusPolynomial(output.a[i], input1.a[i], input2.a[i]);
    }
    subTorusPolynomial(output.b, input1.b, input2.b);
}

// Computes output = scalar - input without requiring a pre-built rlweOne ciphertext.
// Uses subTorus for correct modular reduction under any TORUS_Q.
template<typename TrlweType>
void subTrlweFromConst(TrlweType& output, const TrlweType& input, const Torus scalar) {
    const auto N = output.b.N;
    for (auto j = 0; j < (int)output.a.size(); j++) {
        for (auto i = 0; i < N; i++) {
            output.a[j].coeffs[i] = subTorus(TORUS_Q, Torus(0), input.a[j].coeffs[i]);
        }
    }
    output.b.coeffs[0] = subTorus(TORUS_Q, scalar, input.b.coeffs[0]);
    for (auto i = 1; i < N; i++) {
        output.b.coeffs[i] = subTorus(TORUS_Q, Torus(0), input.b.coeffs[i]);
    }
}

// template<typename TrlweDftType>
// void addTrlweNtt(TrlweDftType& output, const TrlweDftType& input1, const TrlweDftType& input2) {
//     for (auto i = 0; i < output.a.size(); i++) {
//         addNttPolynomial(output.a[i], input1.a[i], input2.a[i]);
//     }
//     addNttPolynomial(output.b, input1.b, input2.b);
// }
template<typename TrlweDftType, typename... TrlweDftArgs>
void addTrlweNtt(TrlweDftType& output, const TrlweDftArgs&... inputs) {
    for (auto i = 0; i < output.a.size(); i++) {
        NttHexl::addNttPolynomial(output.a[i], inputs.a[i]...);
    }
    NttHexl::addNttPolynomial(output.b, inputs.b...);
}

template<typename TrlweDftType>
void subTrlweNtt(TrlweDftType& output, const TrlweDftType& input1, const TrlweDftType& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        subNttPolynomial(output.a[i], input1.a[i], input2.a[i]);
    }
    subNttPolynomial(output.b, input1.b, input2.b);
}

template <typename T>
void resetTrlweToZero(std::vector<T>& a, T& b) {
    std::fill(a.begin(), a.end(), T(b.N, 0));
    std::fill(b.coeffs.begin(), b.coeffs.end(), 0);
}

template<typename TrlweType>
void accumulateTrlwe(TrlweType& accum, const TrlweType& trlwe) {
    for (auto i = 0; i < accum.a.size(); i++) {
        accumulateTorusPolynomial(accum.a[i], trlwe.a[i]);
    }
    accumulateTorusPolynomial(accum.b, trlwe.b);
}

template<typename TrlweType, typename U>
void accumulateTrlweModP(TrlweType& accum, const TrlweType& trlwe, const U p) {
    for (auto i = 0; i < accum.a.size(); i++) {
        accumulatePolynomialModP(accum.a[i], trlwe.a[i], p);
    }
    accumulatePolynomialModP(accum.b, trlwe.b, p);
}

template<typename TrlweTypeA, typename TrlweTypeB>
void decompTrlweCrt(std::vector<TrlweTypeA>& out, const TrlweTypeB& in, const YatfheParameters& param) {
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
void decompTrlweCrtNO(TrlweTypeA& out, const TrlweTypeB& in, const YatfheParameters& param) {
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
void decompTrlweMcrt(std::vector<TrlweTypeA>& out, const TrlweTypeB& in, const YatfheParameters& param) {
    const size_t d = param.d;
    const size_t N = param.N;

    int64_t taoU[NUM_PRIMES];
    int     qd[NUM_PRIMES];
    for (size_t di = 0; di < d; di++) {
        taoU[di] = param.taoU[di];
        qd[di]   = param.qd[di];
    }

    using ACoeffPtr = decltype(out[0].a[0].coeffs.data());
    using BCoeffPtr = decltype(out[0].b.coeffs.data());

    for (size_t ki = 0; ki < static_cast<size_t>(param.k); ki++) {
        const auto* inA = in.a[ki].coeffs.data();
        ACoeffPtr outA[NUM_PRIMES];
        for (size_t di = 0; di < d; di++) outA[di] = out[di].a[ki].coeffs.data();

        for (size_t j = 0; j < N; j++) {
            const int64_t val = inA[j];
            for (size_t di = 0; di < d; di++) {
                outA[di][j] = static_cast<int8_t>(longModP(taoU[di] * val, qd[di]));
            }
        }
    }

    const auto* inB = in.b.coeffs.data();
    BCoeffPtr outB[NUM_PRIMES];
    for (size_t di = 0; di < d; di++) outB[di] = out[di].b.coeffs.data();

    for (size_t j = 0; j < N; j++) {
        const int64_t val = inB[j];
        for (size_t di = 0; di < d; di++) {
            outB[di][j] = static_cast<int8_t>(longModP(taoU[di] * val, qd[di]));
        }
    }
}

template<typename TrlweType>
void trlweMcrtToCrt(std::vector<TrlweType>& trlwe, const YatfheParameters& param) {
    const auto dh = param.dh;
    const size_t dl = param.dl;
    const size_t N  = param.N;

    // Hoist low-prime constants into stack arrays so all dl values are available
    // simultaneously when the loop order is (k, j, d) rather than (d, k, j).
    int64_t taoUInv[NUM_LOW_PRIMES];
    int     ql[NUM_LOW_PRIMES];
    for (size_t d = 0; d < dl; d++) {
        taoUInv[d] = param.taoUInv[dh + d];
        ql[d]      = param.ql[d];
    }

    using CoeffPtr = decltype(trlwe[0].a[0].coeffs.data());

    for (size_t k = 0; k < static_cast<size_t>(param.k); k++) {
        CoeffPtr pA[NUM_LOW_PRIMES];
        for (size_t d = 0; d < dl; d++) pA[d] = trlwe[dh + d].a[k].coeffs.data();

        for (size_t j = 0; j < N; j++) {
            for (size_t d = 0; d < dl; d++) {
                pA[d][j] = static_cast<int8_t>(longModP(taoUInv[d] * pA[d][j], ql[d]));
            }
        }
    }

    CoeffPtr pB[NUM_LOW_PRIMES];
    for (size_t d = 0; d < dl; d++) pB[d] = trlwe[dh + d].b.coeffs.data();

    for (size_t j = 0; j < N; j++) {
        for (size_t d = 0; d < dl; d++) {
            pB[d][j] = static_cast<int8_t>(longModP(taoUInv[d] * pB[d][j], ql[d]));
        }
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void recompTrlweApproxCrt(TrlweTypeA& out, std::vector<TrlweTypeB>& inMCRT, const YatfheParameters& param) {
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
void recompTrlweCrt(TrlweTypeA& out, std::vector<TrlweTypeB>& inCRT, const YatfheParameters& param) {
    const size_t d   = param.d;
    const size_t N   = param.N;
    const auto   qCRT = param.qCRT;

    // Hoist gadget vector to avoid param pointer chasing in the hot loop.
    long z[NUM_PRIMES];
    for (size_t di = 0; di < d; di++) z[di] = param.z[di];

    // Cache one raw input pointer per prime outside the j loop to remove the
    // inCRT[d].a[k].coeffs chain (three pointer dereferences) from the inner body.
    using InCoeffPtr = decltype(inCRT[0].a[0].coeffs.data());

    for (size_t k = 0; k < static_cast<size_t>(param.k); k++) {
        auto* outA = out.a[k].coeffs.data();
        InCoeffPtr inA[NUM_PRIMES];
        for (size_t di = 0; di < d; di++) inA[di] = inCRT[di].a[k].coeffs.data();

        for (size_t j = 0; j < N; j++) {
            long tmp = 0;
            for (size_t di = 0; di < d; di++) tmp += inA[di][j] * z[di];
            outA[j] = static_cast<Torus>(longModP(tmp, qCRT));
        }
    }

    auto* outB = out.b.coeffs.data();
    InCoeffPtr inB[NUM_PRIMES];
    for (size_t di = 0; di < d; di++) inB[di] = inCRT[di].b.coeffs.data();

    for (size_t j = 0; j < N; j++) {
        long tmp = 0;
        for (size_t di = 0; di < d; di++) tmp += inB[di][j] * z[di];
        outB[j] = static_cast<Torus>(longModP(tmp, qCRT));
    }
}

template<typename TrlweTypeA, typename TrlweTypeB>
void recompTrlweCrtNO(TrlweTypeA& out, TrlweTypeB& inCRT, const YatfheParameters& param) {
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

template<typename TrlweType>
void clearTrlwe(TrlweType& obj) {
    for (auto & poly : obj.a) {
        std::fill(poly.coeffs.begin(), poly.coeffs.end(), 0);
    }
    std::fill(obj.b.coeffs.begin(), obj.b.coeffs.end(), 0);
}

// CKKS-style mod-down: divide every coefficient by 2^shift with round-to-nearest.
// Reduces per-ciphertext noise by ~2^shift at the cost of scaling the message down
// by the same factor; callers must decode with torusBase << shift afterwards.
inline void modDownTrlwe(Trlwe& trlwe, const int shift) {
    if (shift <= 0) return;
    const int64_t half = static_cast<int64_t>(1) << (shift - 1);
    for (auto& poly : trlwe.a) {
        for (Torus& c : poly.coeffs) {
            c = static_cast<Torus>(c >> shift);
        }
    }
    for (Torus& c : trlwe.b.coeffs) {
        c = static_cast<Torus>(c >> shift);
    }
}

void genTrlweKey(TrlweKey& key);

void symEncTrlweSingleSample(Trlwe& trlwe, const TrlweKey& key, Torus mu, int pos);

void symEncTrlweSingleSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, Torus mu, int pos);

void symEncTrlweSingleSampleNttSimple(TrlweDft& trlweDft, const TrlweKey& key, Torus mu, int pos);

void symEncTrlweMultiSample(Trlwe& trlwe, const TrlweKey& key, const std::vector<Torus>& mu);

void symEncTrlweMultiSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const std::vector<Torus>& mu);

void symEncTrlweMultiSampleSimple(TrlweDft& trlweDft, const TrlweKey& key, const vector<Torus>& mu);

void symDecTrlweToDouble(DoublePolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweToTorus(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweToInt(IntPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, int torusBase);

void symDecTrlweToIntNtt(IntPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRounding(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key);

void symDecTrlweNtt(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, int torusBase);

void symDecTrlweWoRoundingNtt(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param);

void gadgetDecomposeTrlweA(vector<vector<DecompPolynomial>>& output, const vector<TorusPolynomial>& a, const YatfheParameters& param);
void gadgetDecomposeTrlweB(vector<DecompPolynomial>& output, const TorusPolynomial& b, const YatfheParameters& param);

void recomposeTrlwe(Trlwe& output, const DecomposedTrlwe& input, const YatfheParameters& param);

void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, int index);

void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey);

void rotateTrlwe(Trlwe& res, const Trlwe& input, int a);

void rotateAccumulateTrlwe(Trlwe& accum, const Trlwe& input, int aTrue, int isWrap);

void rotateTrlweNtt(TrlweDft& res, const TrlweDft& input, const int r);

void rotateTrlweMinusOne(Trlwe& res, const Trlwe& input, int a);

void rotateTrlweMinusOneBPlusOne(Trlwe& res, TorusPolynomial& b, const Trlwe& input, int a, Torus one);

void rotateTrlweMinusOneNtt(TrlweDft& res, const TrlweDft& input, const int r);

void rotateTrlweMinusOneBPlusOneNtt(TrlweDft& res, const TrlweDft& input, const int r);

void rotateTrlwe8MinusOne(Trlwe8& res, const Trlwe8& input, int a, int modP);

void copyTrlwe(Trlwe& target, const Trlwe& source, bool copyA, bool copyB);

void rescaleTrlweToNewMod(Trlwe& output, const Trlwe& in, int64_t newMod, int64_t currMod);

void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput);

void multTrlweWithConst(Trlwe& output, const Trlwe& input1, const int scalar);

void multTrlweWithPolyNtt(Trlwe& output, const Trlwe& in, const IntPolynomial& poly, int level, const YatfheParameters& param);

#endif //HLS_YATFHE_TRLWE_H
