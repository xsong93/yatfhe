//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include <random>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

void trlweKeyGen(TrlweKey& key, const int N, const int k) {
    uniform_int_distribution<int> distribution(0, 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < N; j++) {
            key.s[i].coeffs[j] = distribution(rng);
        }
    }
//    printPolyVec(key.s, "TrlweKey");
}

// Trlwe: (X^-b) * (0,...,0,v)
void genNoiselessTrlweSample(Trlwe& accum, const TorusPolynomial& v, const ScaledTlwe& scaledInput) {
    const auto barb = scaledInput.b;
    torusPolynomialRotate(accum.b, -barb, v);

    std::vector<double> t(accum.b.coeffs.size());
    for (int i = 0; i < accum.b.N; i++) {
        t[i] = torus32ToDouble(accum.b.coeffs[i]);
    }
    printArray(t, "b:");
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

// out = (a', b[index])
void extractTlweFromTrlwe(Tlwe& out, const Trlwe& in, const int index) {
    const auto N = in.b.N;
    const auto size = in.a.size();
    for (auto i = 0; i < size; i++) {
        for (auto j = 0; j < N; j++) {
            out.a[i * N + j] = (j <= index) ? (in.a[i].coeffs[index - j]) : (-in.a[i].coeffs[N + index - j]);
        }
    }
    out.b = in.b.coeffs[index];
}

// res = X^a * input - input
void trlweRotateMinusOne(Trlwe& res, const Trlwe& input, const int a) {
    const auto size = input.a.size();
    for (auto i = 0; i < size; i++) {
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

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposition(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto l = param.l;
    const auto radixBits = param.radixBits;
    const auto maskMod = param.maskMod;
    const auto bHalf = param.bHalf;
    const auto offset = genOffset(radixBits, bHalf, l, param.torusBits);
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
        polynomialAddSubOffset(currIn, offset, true);
        for (auto lvl = 0; lvl < l; lvl++) {
            const auto decal = 32 - (lvl + 1) * radixBits;
            for (auto j = 0; j < N; j++) {
                auto& currOut = (row < k) ? output.rlwes[lvl].a[row] : output.rlwes[lvl].b;
                currOut.coeffs[j] = (currIn.coeffs[j] >> decal) & maskMod - bHalf;
            }
        }
        polynomialAddSubOffset(currIn, offset, false);
    }
}
