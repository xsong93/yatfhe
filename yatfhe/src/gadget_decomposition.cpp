//
// Created by Xintong Song on 2024/3/15.
//
#include <iostream>
#include <vector>
#include <numeric>
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/gadget_decomposition.h"

using namespace std;

// offset = B/2 * (2^(torusBits - radixBits) + 2^(torusBits - 2 * radixBits) + ... + 2^(torusBits - l * radixBits))
int genOffset(const int radixBits, const int bHalf, const int l, const int torusBits) {
    int res = 0;
    for (auto i = 1; i <= l; ++i) {
        res += 1 << (torusBits - i * radixBits);
    }
    return res * bHalf;
}

// g = (1/B, ..., 1/B^l), B = 2^radixBits
std::vector<Torus> genGadgetVector(const int radixBits, const int l, const int torusBits) {
    std::vector<Torus> g(l);
    for (auto i = 1; i <= l; i++) {
        g[i - 1] = 1 << (torusBits - i * radixBits); // 1/(B^(i) as Torus: 2^torusBits * 2^(-radixBits*i)
    }
    return g;
}

UnsignedInteger recompose(const std::vector<UnsignedInteger>& digits, const YatfheParameters& param) {
    std::vector<UnsignedInteger> shiftedDigits(digits.size());
    for (auto i = 1; i <= digits.size(); ++i) {
        shiftedDigits[i - 1] = digits[i] << (param.torusBits - i * param.radixBits);
    }
    return accumulate(shiftedDigits.begin(), shiftedDigits.end(), 0u);
}

/**
 * Calculate in * B^-j.
 * @param in The input to decompose.
 * @param param
 * @return
 */
std::vector<Integer> decomposeOverB(const Integer in, const YatfheParameters& param) {
    std::vector<Integer> output(param.ksLevel);
    for (int i = 1; i <= output.size(); ++i) {
        output[i - 1] = in << (param.torusBits - i * param.radixBits);
    }
    return output;
}

/**
 * Calculate signed decomposition g^-1(x).
 * @param res Resulting g^-1(x).
 * @param input The input to decompose.
 * @param param
 */
void signedGadgetDecomposition(vector<Torus>& res, const Torus in, const YatfheParameters& param) {
    vector<Torus> tmp(param.torusBits / param.radixBits);
    auto carry = 0u;
    for (auto i = 0; i < tmp.size(); i++) {
        auto unsignedDigit = ((in >> (i * param.radixBits)) & param.digitMask) + carry;
        auto carryMask = unsignedDigit & param.baseOverTwo;
        auto signedDigit = unsignedDigit - (carryMask << 1);
        carry = carryMask >> (param.radixBase - 1);
        tmp[tmp.size() - i - 1] = signedDigit;
    }
    copy(tmp.begin(), tmp.begin() + param.ksLevel, res.begin());
}

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto l = param.l;
    const auto radixBits = param.radixBits;
    const auto maskMod = param.digitMask;
    const auto baseOverTwo = param.baseOverTwo;
    const auto offset = genOffset(radixBits, baseOverTwo, l, param.torusBits);
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
        polynomialAddSubOffset(currIn, offset, true);
        for (auto lvl = 0; lvl < l; lvl++) {
            const auto decal = param.torusBits - (lvl + 1) * radixBits;
            for (auto j = 0; j < N; j++) {
                auto& currOut = (row < k) ? output.rlwes[lvl].a[row] : output.rlwes[lvl].b;
                currOut.coeffs[j] = (currIn.coeffs[j] >> decal) & maskMod - baseOverTwo;
            }
        }
        polynomialAddSubOffset(currIn, offset, false);
    }
}