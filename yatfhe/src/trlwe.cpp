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
#include "yautil/control_helper.h"

using namespace std;
using namespace NttHexl;

void initTrlweSingleSample(Trlwe& trlwe, const Torus mu, double sigma) {
    initCoeffsWithGaussianNoiseSingleSample(trlwe.b.coeffs, mu, sigma, TORUS_Q);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs, TORUS_MIN, TORUS_MAX);
    }
}

void initTrlweSingleSampleFixedNoise(Trlwe& trlwe, const Torus mu, const Torus noise) {
    for (auto i = 0 ; i < trlwe.N; i++) {
        trlwe.b.coeffs[i] = longModP(mu + noise, TORUS_Q);
    }
    for (auto i = 0 ; i < trlwe.k; i++) {
        for (auto j = 0 ; j < trlwe.N; j++) {
            trlwe.a[i].coeffs[j] = 0;
        }
    }
}

void initTrlweMultiSample(Trlwe& trlwe, const vector<Torus>& mu, double sigma) {
    initCoeffsWithGaussianNoiseMultiSample(trlwe.b.coeffs, mu, sigma, TORUS_Q);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs, TORUS_MIN, TORUS_MAX);
    }
}

void initTrlweMultiSampleFixedNoise(Trlwe& trlwe, const vector<Torus>& mu, const Torus noise) {
    for (auto i = 0 ; i < trlwe.N; i++) {
        trlwe.b.coeffs[i] = longModP(mu[i] + noise, TORUS_Q);
    }
    for (auto i = 0 ; i < trlwe.k; i++) {
        for (auto j = 0 ; j < trlwe.N; j++) {
            trlwe.a[i].coeffs[j] = 0;
        }
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

void symEncTrlweSingleSample(Trlwe& trlwe, const TrlweKey& key, const Torus mu) {
    initTrlweSingleSample(trlwe, mu, key.sigma);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweSingleSampleFixedNoise(Trlwe& trlwe, const TrlweKey& key, const Torus mu, const Torus noise) {
    initTrlweSingleSampleFixedNoise(trlwe, mu, noise);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweMultiSample(Trlwe& trlwe, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSample(trlwe, mu, key.sigma);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweMultiSampleFixedNoise(Trlwe& trlwe, const TrlweKey& key, const vector<Torus>& mu, const Torus noise) {
    initTrlweMultiSampleFixedNoise(trlwe, mu, noise);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweSingleSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const Torus mu) {
    initTrlweSingleSample(trlwe, mu, key.sigma);
    symEncTrlweNtt(trlwe, trlweDft, key);
}

void symEncTrlweMultiSampleNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSample(trlwe, mu, key.sigma);
    symEncTrlweNtt(trlwe, trlweDft, key);
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
        output.coeffs[i] = roundTorusError(output.coeffs[i], torusBase);
    }
}

void symDecTrlweToInt(IntPolynomial& output, const Trlwe& trlwe, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    symDecTrlweWoRounding(tmp, trlwe, key);
    for (auto i = 0 ; i < tmp.N; i++) {
        tmp.coeffs[i] = roundTorusError(tmp.coeffs[i], torusBase);
        output.coeffs[i] = modSwitchFromTorus32(tmp.coeffs[i], torusBase);
    }
}

void symDecTrlweToIntNtt(IntPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    symDecTrlweWoRoundingNtt(tmp, trlweDft, key);
    for (auto i = 0; i < tmp.N; i++) {
        tmp.coeffs[i] = roundTorusError(tmp.coeffs[i], torusBase);
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

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposeTrlwe(DecomposedTrlwe& output, const Trlwe& input, const YatfheParameters& param) {
    const auto k = input.k;
    const auto N = input.b.coeffs.size();
    const auto l = output.l;
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
        for (auto j = 0; j < N; j++) {
            DecomposedData d {l};
//            signedGadgetDecomposition(d, currIn.coeffs[j], param);
            gadgetDecompose(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < l; lvl++) {
                auto& currOut = (row < k) ? output.trlwes[lvl].a[row] : output.trlwes[lvl].b;
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }
}

void gadgetDecomposeTrlweA(vector<vector<DecompPolynomial>>& output, const vector<TorusPolynomial>& a, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto l = param.l;
    for (auto row = 0; row < k; row++) {
        auto& currIn = a[row];
        for (auto j = 0; j < N; j++) {
            DecomposedData d{l};
            gadgetDecompose(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < l; lvl++) {
                auto& currOut = output[lvl][row];
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }
}

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposeTrlweNtt(DecomposedTrlweDft& output, const TrlweDft& input, const YatfheParameters& param) {
    const auto k = input.k;
    const auto N = input.b.coeffs.size();
    const auto l = output.l;
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
        for (auto j = 0; j < N; j++) {
            DecomposedDataDft d {l};
//            signedGadgetDecompositionNtt(d, currIn.coeffs[j], param);
            gadgetDecomposeNtt(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < l; lvl++) {
                auto& currOut = (row < k) ? output.rlweDfts[lvl].a[row] : output.rlweDfts[lvl].b;
                currOut.coeffs[j] = d.value[lvl];
            }
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

// Combine l decomposed TrlweDft a & b into one.
void recomposeTrlweNtt(TrlweDft& output, const DecomposedTrlweDft& input, const YatfheParameters& param) {
    const auto k = output.k;
    const auto N = output.b.coeffs.size();
    const auto l = input.l;
    resetTrlweToZero(output.a, output.b);
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            auto& currIn = (row < k) ? input.rlweDfts[lvl].a[row] : input.rlweDfts[lvl].b;
            auto& currOut = (row < k) ? output.a[row] : output.b;
            for (auto j = 0; j < N; j++) {
                currOut.coeffs[j] = AddUIntMod(currOut.coeffs[j], currIn.coeffs[j] << (param.dftBits - (lvl + 1) * param.radixBits), param.qNtt);
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
    for (auto i = 0; i < input.a.size(); i++) {
        rotateTorusPolynomial(res.a[i], a, input.a[i]);
    }
    rotateTorusPolynomial(res.b, a, input.b);
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
