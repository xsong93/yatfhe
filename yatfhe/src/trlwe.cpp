//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/multi_threading.h"
#include "yautil/tool.h"

using namespace std;
using namespace NttHexl;

namespace {
    void initTrlweSingleSample(Trlwe& trlwe, const Torus mu, const int pos, double sigma) {
        initCoeffsWithTUniformNoiseSingleSample(trlwe.b.coeffs, mu, pos, sigma, TORUS_Q);
        for (auto i = 0 ; i < trlwe.k; i++) {
            initCoeffsViaUniformDistribution(trlwe.a[i].coeffs, TORUS_MIN, TORUS_MAX);
        }
    }

    void initTrlweSingleSampleSimple(TrlweDft& trlweDft, const Torus mu, const int pos, double sigma) {
        TorusPolynomial bTmp{trlweDft.b.N};
        initCoeffsWithTUniformNoiseSingleSample(bTmp.coeffs, mu, pos, sigma, TORUS_Q);
        applyNtt(trlweDft.b, bTmp);
        for (auto i = 0 ; i < trlweDft.k; i++) {
            initNttCoeffsViaUniformDistribution(trlweDft.a[i].coeffs, NTT_MIN, NTT_MAX);
        }
    }

    void initTrlweMultiSample(Trlwe& trlwe, const vector<Torus>& mu, double sigma) {
       initCoeffsWithTUniformNoiseMultiSample(trlwe.b.coeffs, mu, sigma, TORUS_Q);
       for (auto i = 0 ; i < trlwe.k; i++) {
           initCoeffsViaUniformDistribution(trlwe.a[i].coeffs, TORUS_MIN, TORUS_MAX);
       }
    }

    void initTrlweMultiSampleSimple(TrlweDft& trlweDft, const vector<Torus>& mu, double sigma) {
        TorusPolynomial bTmp{trlweDft.b.N};
        initCoeffsWithTUniformNoiseMultiSample(bTmp.coeffs, mu, sigma, TORUS_Q);
        applyNtt(trlweDft.b, bTmp);
        for (auto i = 0 ; i < trlweDft.k; i++) {
            initNttCoeffsViaUniformDistribution(trlweDft.a[i].coeffs, NTT_MIN, NTT_MAX);
        }
    }

    void symEncTrlwe(Trlwe& trlwe, const TrlweKey& key) {
        for (auto i = 0; i < trlwe.k; i++) {
            multTorusPolynomialAcc(trlwe.b, trlwe.a[i], key.s[i]);
        }
    }

    void symEncTrlweNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key) {
        applyNttForAB(trlweDft, trlwe);
        calModularInnerProductNtt(trlweDft.b, trlweDft.a, key.sDft);
        applyIntt(trlwe.b, trlweDft.b);
    }

    void symEncTrlweNttSimple(TrlweDft& trlweDft, const TrlweKey& key) {
        calModularInnerProductNtt(trlweDft.b, trlweDft.a, key.sDft);
    }
}


void genTrlweKey(TrlweKey& key) {
    for (int i = 0; i < key.k; i++) {
        for (int j = 0; j < key.N; j++) {
#ifdef TERNARY
            key.s[i].coeffs[j] = ternaryDistrib(rng);
#else
            key.s[i].coeffs[j] = binaryDistrib(rng);
#endif
        }
        applyNtt(key.sDft[i], key.s[i]);
    }
}

void symEncTrlweSingleSample(Trlwe& trlwe, const TrlweKey& key, const Torus mu, const int pos) {
    initTrlweSingleSample(trlwe, mu, pos, key.errorB);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweSingleSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const Torus mu, const int pos) {
    initTrlweSingleSample(trlwe, mu, pos, key.errorB);
    symEncTrlweNtt(trlwe, trlweDft, key);
}

void symEncTrlweSingleSampleNttSimple(TrlweDft& trlweDft, const TrlweKey& key, const Torus mu, const int pos) {
    initTrlweSingleSampleSimple(trlweDft, mu, pos, key.errorB);
    symEncTrlweNttSimple(trlweDft, key);
}

void symEncTrlweMultiSample(Trlwe& trlwe, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSample(trlwe, mu, key.errorB);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweMultiSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSample(trlwe, mu, key.errorB);
    symEncTrlweNtt(trlwe, trlweDft, key);
}

void symEncTrlweMultiSampleSimple(TrlweDft& trlweDft, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSampleSimple(trlweDft, mu, key.errorB);
    symEncTrlweNttSimple(trlweDft, key);
}

void symDecTrlweToDouble(DoublePolynomial& output, const Trlwe& trlwe, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    TorusPolynomial innerProduct {output.N};
    for (auto i = 0; i < trlwe.k; i++) {
        multTorusPolynomialAcc(innerProduct, trlwe.a[i], key.s[i]);
    }
    subTorusPolynomial(tmp, trlwe.b, innerProduct);
    torusPolyToDoublePoly(output, tmp);
    roundErrorDoublePoly(output, torusBase);
}

void symDecTrlweToTorus(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, const int torusBase) {
    symDecTrlweWoRounding(output, trlwe, key);
    for (auto i = 0 ; i < output.N; i++) {
        output.coeffs[i] = roundTorus32Error(output.coeffs[i], torusBase);
    }
}

void symDecTrlweToInt(IntPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    symDecTrlweWoRounding(tmp, trlwe, key);
    for (auto i = 0 ; i < tmp.N; i++) {
        tmp.coeffs[i] = roundTorus32Error(tmp.coeffs[i], torusBase);
        output.coeffs[i] = modSwitchFromTorus32(tmp.coeffs[i], torusBase);
    }
}

void symDecTrlweToIntNtt(IntPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    symDecTrlweWoRoundingNtt(tmp, trlweDft, key);
    for (auto i = 0; i < tmp.N; i++) {
        tmp.coeffs[i] = roundTorus32Error(tmp.coeffs[i], torusBase);
        output.coeffs[i] = modSwitchFromTorus32(tmp.coeffs[i], torusBase);
    }
}

void symDecTrlweWoRounding(TorusPolynomial& output, const Trlwe& trlwe, const TrlweKey& key) {
    for (auto i = 0; i < trlwe.k; i++) {
        multTorusPolynomialAcc(output, trlwe.a[i], key.s[i]);
    }
    subTorusPolynomial(output, trlwe.b, output);
}

void symDecTrlweNtt(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    NttPolynomial innerProduct {trlweDft.b.N};
    NttPolynomial res {trlweDft.b.N};
    calModularInnerProductNtt(innerProduct, trlweDft.a, key.sDft);
    subNttPolynomial(res, trlweDft.b, innerProduct);
    applyIntt(tmp, res);
    torusPolyToDoublePoly(output, tmp);
    roundErrorDoublePoly(output, torusBase);
}

void symDecTrlweWoRoundingNtt(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key) {
    NttPolynomial tmp {trlweDft.b.N};
    calModularInnerProductNtt(tmp, trlweDft.a, key.sDft);
    subNttPolynomial(tmp, trlweDft.b, tmp);
    applyIntt(output, tmp);
}

// Trlwe: (X^-b) * (0,...,0,v)
void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput) {
    const auto barb = scaledInput.b;
    rotateTorusPolynomial(accum.b, -barb, v);
}

namespace {
    // Straight-line specialisation of signedGadgetDecomposition for a
    // compile-time digit count L. The generic version leaves the coefficient loop
    // with two consecutive inner loops (digit extraction, then write-out), which
    // GCC refuses to vectorise,
    //
    // The carry is per-coefficient (reset every j), so coefficients stay
    // independent and vectorising across j is exact.
    template<int L>
    void decomposeRowUnrolled(Torus* const* outPtr, const Torus* in, const int N,
                              const int radixBits, const int torusBits) {
        const Torus B = static_cast<Torus>(1) << radixBits;
        const Torus halfB = B >> 1;
        const int shift = torusBits - L * radixBits;
        const UnsignedInteger round =
            (shift > 0) ? (static_cast<UnsignedInteger>(1) << (shift - 1)) : 0;

        Torus* out[L];
        for (int lvl = 0; lvl < L; lvl++) out[lvl] = outPtr[lvl];

        for (int j = 0; j < N; j++) {
            const UnsignedInteger u = static_cast<UnsignedInteger>(in[j]) + round;
            Torus carry = 0;
            for (int lvl = L - 1; lvl >= 0; --lvl) {   // constant trip count -> fully unrolled
                const UnsignedInteger window =
                    (u >> (torusBits - (lvl + 1) * radixBits)) & static_cast<UnsignedInteger>(B - 1);
                Torus digit = static_cast<Torus>(window) + carry;
                carry = (digit >= halfB);
                digit -= carry * B;
                out[lvl][j] = digit;                   // sign is always +1, so the
            }                                          // generic "* d.sign" is a no-op
        }
    }
}

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param) {
    const auto k = input.k;
    const int N = static_cast<int>(input.b.coeffs.size());
    const auto l = output.l;
    DecomposedData d{l};
    // Per-row output pointers, resolved once instead of walking
    // output.trlwes[lvl].a[row] again for every coefficient.
    std::vector<Torus*> outPtr(l);
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
        for (auto lvl = 0; lvl < l; lvl++) {
            auto& currOut = (row < k) ? output.trlwes[lvl].a[row] : output.trlwes[lvl].b;
            outPtr[lvl] = currOut.coeffs.data();
        }
        switch (l) {
            case 1: decomposeRowUnrolled<1>(outPtr.data(), currIn.coeffs.data(), N, param.radixBits, param.torusBits); break;
            case 2: decomposeRowUnrolled<2>(outPtr.data(), currIn.coeffs.data(), N, param.radixBits, param.torusBits); break;
            case 3: decomposeRowUnrolled<3>(outPtr.data(), currIn.coeffs.data(), N, param.radixBits, param.torusBits); break;
            case 4: decomposeRowUnrolled<4>(outPtr.data(), currIn.coeffs.data(), N, param.radixBits, param.torusBits); break;
            default:
                for (auto j = 0; j < N; j++) {
                    signedGadgetDecomposition(d, currIn.coeffs[j], param);
                    for (auto lvl = 0; lvl < l; lvl++) {
                        outPtr[lvl][j] = d.value[lvl] * d.sign;
                    }
                }
        }
    }
}

void gadgetDecomposeTrlweA(vector<vector<DecompPolynomial>>& output, const vector<TorusPolynomial>& a, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto l = param.l;
    DecomposedData d{l};
    for (auto row = 0; row < k; row++) {
        auto& currIn = a[row];
        for (auto j = 0; j < N; j++) {
            signedGadgetDecomposition(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < l; lvl++) {
                auto& currOut = output[lvl][row];
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }
}

void gadgetDecomposeTrlweB(vector<DecompPolynomial>& output, const TorusPolynomial& b, const YatfheParameters& param) {
    const auto N = param.N;
    const auto l = param.l;
    DecomposedData d{l};
    for (auto j = 0; j < N; j++) {
        signedGadgetDecomposition(d, b.coeffs[j], param);
        for (auto lvl = 0; lvl < l; lvl++) {
            output[lvl].coeffs[j] = d.value[lvl] * d.sign;
        }
    }
}

// Combine l decomposed Trlwe a & b into one.
void recomposeTrlwe(Trlwe& output, const DecomposedTrlwe& input, const YatfheParameters& param) {
    const auto k = output.k;
    const auto N = output.b.coeffs.size();
    const auto l = input.l;
    resetTrlweToZero(output.a, output.b);
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            auto& currIn = (row < k) ? input.trlwes[lvl].a[row] : input.trlwes[lvl].b;
            auto& currOut = (row < k) ? output.a[row] : output.b;
            for (auto j = 0; j < N; j++) {
                currOut.coeffs[j] += currIn.coeffs[j] << (param.torusBits - (lvl + 1) * param.radixBits);
            }
        }
    }
}

// out = (a', b[index])
void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, const int index) {
//    tlweCLear(out);
    const auto N = in.b.N;
    for (auto i = 0; i < in.k; i++) {
        for (auto j = 0; j < N; j++) {
            out.a[i * N + j] = modSwitchFromTorusGeneral(j <= index ? in.a[i].coeffs[index - j] : -in.a[i].coeffs[N + index - j], LWE_Q, TORUS_Q);
        }
    }
    out.b = modSwitchFromTorusGeneral(in.b.coeffs[index], LWE_Q, TORUS_Q);
}

// flatten a trlwe key as a tlwe key
void convertTrlweKeyToTlweKey(TlweKey& tlweKey, const TrlweKey& trlweKey) {
    for (auto i = 0; i < trlweKey.k; i++) {
        const auto N = trlweKey.s[i].N;
        for (auto j = 0; j < N; j++) {
            tlweKey.s[i * N + j] = trlweKey.s[i].coeffs[j];
        }
    }
}

void rotateTrlwe(Trlwe& res, const Trlwe& input, const int a) {
    if (a % (input.b.N * 2) == 0) return;
    for (auto i = 0; i < input.a.size(); i++) {
        rotateTorusPolynomial(res.a[i], a, input.a[i]);
    }
    rotateTorusPolynomial(res.b, a, input.b);
}

void rotateAccumulateTrlwe(Trlwe& accum, const Trlwe& input, const int aTrue, const int isWrap) {
    for (auto i = 0; i < static_cast<int>(input.a.size()); i++) {
        rotateAccumulateTorusPolynomial(accum.a[i], aTrue, isWrap, input.a[i]);
    }
    rotateAccumulateTorusPolynomial(accum.b, aTrue, isWrap, input.b);
}

void rotateTrlweNtt(TrlweDft& res, const TrlweDft& input, const int r) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateNttPolynomial(res.a[i], input.a[i], r);
    }
    rotateNttPolynomial(res.b, input.b, r);
}

// res = X^a * input - input
void rotateTrlweMinusOne(Trlwe& res, const Trlwe& input, const int a) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateTorusPolynomialMinusOne(res.a[i], a, input.a[i]);
    }
    rotateTorusPolynomialMinusOne(res.b, a, input.b);
}

void rotateTrlweMinusOneBPlusOne(Trlwe& res, TorusPolynomial& b, const Trlwe& input, const int a, const Torus one) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateTorusPolynomialMinusOne(res.a[i], a, input.a[i]);
    }
    rotateTorusPolynomialMinusOne(b, a, input.b);
    b.coeffs[0] = addTorus(TORUS_Q, b.coeffs[0], one);
}

void rotateTrlweMinusOneNtt(TrlweDft& res, const TrlweDft& input, const int r) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateNttPolynomialMinusOne(res.a[i], input.a[i], r);
    }
    rotateNttPolynomialMinusOne(res.b, input.b, r);
}

void rotateTrlweMinusOneBPlusOneNtt(TrlweDft& res, const TrlweDft& input, const int r) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateNttPolynomialMinusOne(res.a[i], input.a[i], r);
    }
    rotateNttPolynomialMinusOne(res.b, input.b, r);
    const auto q = NttHexl::getNttHexl().GetModulus();
    for (auto& coeff : res.b.coeffs) {
        coeff += 1;
        if (coeff >= q) coeff -= q;
    }
}

void rotateTrlwe8MinusOne(Trlwe8& res, const Trlwe8& input, const int a, int modP) {
    for (auto i = 0; i < input.a.size(); i++) {
        rotateInt8PolynomialMinusOne(res.a[i], a, input.a[i], modP);
    }
    rotateInt8PolynomialMinusOne(res.b, a, input.b, modP);
}

void copyTrlwe(Trlwe& target, const Trlwe& source, const bool copyA, const bool copyB) {
    const auto size = source.a.size();
    const auto N = source.b.N;
    for (auto j = 0; j < N; j++) {
        if (copyA) {
            for (auto i = 0; i < size; i++) {
                target.a[i].coeffs[j] = source.a[i].coeffs[j];
            }
        }
        if (copyB) {
            target.b.coeffs[j] = source.b.coeffs[j];
        }
    }
}

void rescaleTrlweToNewMod(Trlwe& output, const Trlwe& in, const int64_t newMod, const int64_t currMod) {
    for (auto i = 0; i < in.k; i++) {
        auto& aIn = in.a[i];
        auto& aOut = output.a[i];
        for (auto j = 0; j < aIn.N; j++) {
            aOut.coeffs[j] = modSwitchFromTorusGeneral(aIn.coeffs[j], newMod, currMod);
        }
    }
    for (auto j = 0; j < in.b.N; j++) {
        output.b.coeffs[j] = modSwitchFromTorusGeneral(in.b.coeffs[j], newMod, currMod);
    }
}

void multTrlweWithConst(Trlwe& output, const Trlwe& input1, const int scalar) {
    if (scalar == 0) {
        clearTrlwe(output);
        return;
    }

    // // decompose the scalar on base 2
    // std::vector<int> mults;
    // auto cop = scalar;
    // while (cop >= 2) {
    //     mults.emplace_back(2);
    //     cop = cop % 2;
    // }
    // if (cop == 1) {
    //     mults.emplace_back(1);
    // }

    // decompose the scalar on base 1
    for (auto i = 0; i < scalar; i++) {
        accumulateTrlwe(output, input1);
    }
    //
    // Trlwe tmp{input1.k, input1.N};
    // for (const auto m : mults) {
    //     for (auto i = 0; i < output.k; i++) {
    //         for (int j = 0; j < output.N; j++) {
    //             tmp.a[i].coeffs[j] = multTorus(TORUS_Q, input1.a[i].coeffs[j], m);
    //         }
    //     }
    //     for (int j = 0; j < output.N; j++) {
    //         tmp.b.coeffs[j] = multTorus(TORUS_Q, input1.b.coeffs[j], m);
    //     }
    //     accumulateTrlwe(output, tmp);
    // }
}

// no GD, direct NTT mult
void multTrlweWithPolyNtt(Trlwe& output, const Trlwe& in, const IntPolynomial& poly, const YatfheParameters& param) {
    const auto K = param.k;
    TrlweDft inDft{param};
    TrlweDft resDft{param};
    NttPolynomial polyDft{param.N};

    // ntt
    applyNttForAB(inDft, in);
    applyNtt(polyDft, poly);

    // mult, resDft starts at zero so the accumulating product yields in * poly
    for (auto k = 0; k < K; k++) {
        calModularInnerProductNtt(resDft.a[k], inDft.a[k], polyDft);
    }
    calModularInnerProductNtt(resDft.b, inDft.b, polyDft);

    //intt
    applyInttForAB(output, resDft);
}

// Worst-case magnitude of the extra additive error multTrlweWithPolyNtt()
// introduces for this plaintext, in torus units.
int64_t directNttWrapNoise(const IntPolynomial& poly, const YatfheParameters& param) {
    __int128 norm = 0;
    for (const auto c : poly.coeffs) {
        norm += c < 0 ? -static_cast<int64_t>(c) : static_cast<int64_t>(c);
    }
    // centred residue of qNtt mod TORUS_Q: what one wrap costs
    auto delta = static_cast<int64_t>(param.qNtt % static_cast<uint64_t>(TORUS_Q));
    if (delta > TORUS_Q / 2) {
        delta -= TORUS_Q;
    }
    if (delta < 0) {
        delta = -delta;
    }
    // |wraps| <= |true coeff| / qNtt + 1/2, and |true coeff| <= (TORUS_Q / 2) * ||poly||_1
    const __int128 wraps = (static_cast<__int128>(TORUS_Q / 2) * norm) / static_cast<int64_t>(param.qNtt) + 1;
    const __int128 noise = wraps * delta;
    return noise > INT64_MAX ? INT64_MAX : static_cast<int64_t>(noise);
}

// Whether the wrap noise alone stays inside the rounding margin of one message
// step. This is the full margin -- a caller whose ciphertext already carries
// significant noise should compare directNttWrapNoise() against its own budget.
bool isPolyDirectNttSafe(const IntPolynomial& poly, const YatfheParameters& param) {
    return directNttWrapNoise(poly, param) < TORUS_Q / (2 * param.torusBase);
}
