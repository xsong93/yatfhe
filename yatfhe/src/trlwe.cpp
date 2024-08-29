//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"

using namespace std;

void initTrlweSingleSample(Trlwe& trlwe, const Torus mu, double sigma) {
    initCoeffsWithGaussianNoiseSingleSample(trlwe.b.coeffs, mu, sigma);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs);
    }
}

void initTrlweMultiSample(Trlwe& trlwe, const vector<Torus>& mu, double sigma) {
    initCoeffsWithGaussianNoiseMultiSample(trlwe.b.coeffs, mu, sigma);
    for (auto i = 0 ; i < trlwe.k; i++) {
        initCoeffsViaUniformDistribution(trlwe.a[i].coeffs);
    }
}

void symEncTrlwe(Trlwe& trlwe, const TrlweKey& key) {
    for (auto i = 0; i < trlwe.k; i++) {
        polynomialMulAccNaive(trlwe.b, trlwe.a[i], key.s[i]);
    }
}

void symEncTrlweNtt(Trlwe& trlwe, TrlweDft& trlweDft, const TrlweKey& key) {
    applyNttForAB(trlweDft, trlwe);
    calModularInnerProductNtt(trlweDft.b, trlweDft.a, key.sDft);
    applyIntt(trlwe.b, trlweDft.b);
}

void trlweKeyGen(TrlweKey& key) {
    for (int i = 0; i < key.k; i++) {
        for (int j = 0; j < key.N; j++) {
            key.s[i].coeffs[j] = binaryDistrib(rng);
        }
        applyNtt(key.sDft[i], key.s[i]);
    }
//    printPolyVec(key.s, "TrlweKey");
}

void symEncTrlweSingleSample(Trlwe& trlwe, const TrlweKey& key, const Torus mu) {
    initTrlweSingleSample(trlwe, mu, key.sigma);
    symEncTrlwe(trlwe, key);
}

void symEncTrlweMultiSample(Trlwe& trlwe, const TrlweKey& key, const vector<Torus>& mu) {
    initTrlweMultiSample(trlwe, mu, key.sigma);
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
        polynomialMulAccNaive(innerProduct, trlwe.a[i], key.s[i]);
    }
    polynomialSub(tmp, trlwe.b, innerProduct);
    torusPolyToDoublePoly(output, tmp);
    roundErrorPoly(output, torusBase);
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
        polynomialMulAccNaive(output, trlwe.a[i], key.s[i]);
    }
    polynomialSub(output, trlwe.b, output);
}

void symDecTrlweNtt(DoublePolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key, const int torusBase) {
    TorusPolynomial tmp {output.N};
    LagrangePolynomial innerProduct {trlweDft.b.N};
    LagrangePolynomial res {trlweDft.b.N};
    calModularInnerProductNtt(innerProduct, trlweDft.a, key.sDft);
    lagrangePolynomialSub(res, trlweDft.b, innerProduct);
    applyIntt(tmp, res);
    torusPolyToDoublePoly(output, tmp);
    roundErrorPoly(output, torusBase);
}

void symDecTrlweWoRoundingNtt(TorusPolynomial& output, const TrlweDft& trlweDft, const TrlweKey& key) {
    LagrangePolynomial tmp {trlweDft.b.N};
    calModularInnerProductNtt(tmp, trlweDft.a, key.sDft);
    lagrangePolynomialSub(tmp, trlweDft.b, tmp);
    applyIntt(output, tmp);
}

// Trlwe: (X^-b) * (0,...,0,v)
void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput) {
    const auto barb = scaledInput.b;
    torusPolynomialRotate(accum.b, -barb, v);
}

void trlweAdd(Trlwe& output, const Trlwe& input1, const Trlwe& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        polynomialAdd(output.a[i], input1.a[i], input2.a[i]);
    }
    polynomialAdd(output.b, input1.b, input2.b);
}

void trlweSub(Trlwe& output, const Trlwe& input1, const Trlwe& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        polynomialSub(output.a[i], input1.a[i], input2.a[i]);
    }
    polynomialSub(output.b, input1.b, input2.b);
}

void trlweAddNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        lagrangePolynomialAdd(output.a[i], input1.a[i], input2.a[i]);
    }
    lagrangePolynomialAdd(output.b, input1.b, input2.b);
}

void trlweSubNtt(TrlweDft& output, const TrlweDft& input1, const TrlweDft& input2) {
    for (auto i = 0; i < output.a.size(); i++) {
        lagrangePolynomialSub(output.a[i], input1.a[i], input2.a[i]);
    }
    lagrangePolynomialSub(output.b, input1.b, input2.b);
}

/**
 * accum.a += tlwe.a, accum.b += tlwe.b
 * */
void trlweAccumulate(Trlwe& accum, const Trlwe& tlwe) {
    for (auto i = 0; i < accum.a.size(); i++) {
        polynomialAccumulate(accum.a[i], tlwe.a[i]);
    }
    polynomialAccumulate(accum.b, tlwe.b);
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
                auto& currOut = (row < k) ? output.rlwes[lvl].a[row] : output.rlwes[lvl].b;
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
    trlweSetZero(output.a, output.b);
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            auto& currIn = (row < k) ? input.rlwes[lvl].a[row] : input.rlwes[lvl].b;
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
    trlweSetZero(output.a, output.b);
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            auto& currIn = (row < k) ? input.rlweDfts[lvl].a[row] : input.rlweDfts[lvl].b;
            auto& currOut = (row < k) ? output.a[row] : output.b;
            for (auto j = 0; j < N; j++) {
                currOut.coeffs[j] = modAdd(currOut.coeffs[j], currIn.coeffs[j] << (param.dftBits - (lvl + 1) * param.radixBits));
            }
        }
    }
}

// out = (a', b[index])
void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, const int index) {
    tlweCLear(out);
    const auto N = in.b.N;
    for (auto i = 0; i < in.k; i++) {
        for (auto j = 0; j < N; j++) {
            out.a[i * N + j] = (j <= index) ? (in.a[i].coeffs[index - j]) : (-in.a[i].coeffs[N + index - j]);
        }
    }
    out.b = in.b.coeffs[index];
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

void trlweRotate(Trlwe& res, const Trlwe& input, const int a) {
    for (auto i = 0; i < input.a.size(); i++) {
        torusPolynomialRotate(res.a[i], a, input.a[i]);
    }
    torusPolynomialRotate(res.b, a, input.b);
}

// res = X^a * input - input
void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, const int a) {
    for (auto i = 0; i < input.a.size(); i++) {
        torusPolynomialRotateMinusOne(res.a[i], a, input.a[i]);
    }
    torusPolynomialRotateMinusOne(res.b, a, input.b);
}

// todo: torusPolynomialRotateMinusOne mod q
void trlweRotateMinusOneModQ(Trlwe& res, const Trlwe& input, const int a, const int q) {
    for (auto i = 0; i < input.a.size(); i++) {
        torusPolynomialRotateMinusOne(res.a[i], a, input.a[i]);
    }
    torusPolynomialRotateMinusOne(res.b, a, input.b);
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
