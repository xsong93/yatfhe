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
#include "yautil/tool.h"

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

void gadgetDecompose(DecomposedData& out, const Integer in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    UnsignedInteger tmp = (out.sign == 1) ? in : -in;
    UnsignedInteger mask = ((1 << param.radixBits) - 1) << (param.torusBits - param.radixBits);
    for (auto i = 1; i <= param.ksLevel; i++) {
        out.value[i - 1] = (mask & tmp) >> (param.torusBits - i * param.radixBits);
        mask >>= param.radixBits;
    }
}

Integer recompose(const DecomposedData& digits, const YatfheParameters& param) {
    Integer res {0};
    for (auto i = 1; i <= digits.value.size(); ++i) {
        res += digits.value[i - 1] << (param.torusBits - i * param.radixBits);
    }
    return res * digits.sign;
}

/**
 * Calculate in * B^-j.
 * @param in The input to decompose.
 * @param param
 * @return
 */
std::vector<Integer> decomposeOverB(const Binary in, const YatfheParameters& param) {
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
void signedGadgetDecomposition(DecomposedData& out, const Integer in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    UnsignedInteger unsignedIn = (out.sign == 1) ? in : -in;
    vector<UnsignedInteger> tmp(param.torusBits / param.radixBits);
    UnsignedInteger carry = 0;
    for (auto i = 0; i < tmp.size(); i++) {
        auto unsignedDigit = ((unsignedIn >> (i * param.radixBits)) & param.digitMask) + carry;
        auto carryMask = unsignedDigit & param.baseOverTwo;
        auto signedDigit = unsignedDigit - (carryMask << 1);
        carry = carryMask >> (param.radixBits - 1);
        tmp[tmp.size() - i - 1] = signedDigit;
    }
    copy(tmp.begin(), tmp.begin() + param.ksLevel, out.value.begin());
}

// G^-1 * Trlwe = DecomposedTrlwe
void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto l = param.l;
    const auto radixBits = param.radixBits;
    const auto maskMod = param.digitMask;
    const auto baseOverTwo = param.baseOverTwo;
//    const auto offset = genOffset(radixBits, baseOverTwo, l, param.torusBits);
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? input.a[row] : input.b;
//        polynomialAddSubOffset(currIn, offset, true);
        for (auto lvl = 0; lvl < l; lvl++) {
//            const auto decal = param.torusBits - (lvl + 1) * radixBits;
            for (auto j = 0; j < N; j++) {
                auto& currOut = (row < k) ? output.rlwes[lvl].a[row] : output.rlwes[lvl].b;
//                currOut.coeffs[j] = (currIn.coeffs[j] >> decal) & maskMod - baseOverTwo;
//                signedGadgetDecomposition(, currIn.coeffs[j], param);
            }
        }
//        polynomialAddSubOffset(currIn, offset, false);
    }
}

//// G^-1 * Trlwe = DecomposedTrlwe
//void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param) {
//    const auto k = param.k;
//    const auto N = param.N;
//    const auto l = param.l;
//    const auto radixBits = param.radixBits;
//    const auto maskMod = param.digitMask;
//    const auto baseOverTwo = param.baseOverTwo;
//    const auto offset = genOffset(radixBits, baseOverTwo, l, param.torusBits);
//    for (auto row = 0; row < k + 1; row++) {
//        auto& currIn = (row < k) ? input.a[row] : input.b;
//        polynomialAddSubOffset(currIn, offset, true);
//        for (auto lvl = 0; lvl < l; lvl++) {
//            const auto decal = param.torusBits - (lvl + 1) * radixBits;
//            for (auto j = 0; j < N; j++) {
//                auto& currOut = (row < k) ? output.rlwes[lvl].a[row] : output.rlwes[lvl].b;
//                currOut.coeffs[j] = (currIn.coeffs[j] >> decal) & maskMod - baseOverTwo;
//            }
//        }
//        polynomialAddSubOffset(currIn, offset, false);
//    }
//}