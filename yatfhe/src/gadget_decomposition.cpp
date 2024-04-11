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
#include "yatfhe/ntt.h"
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
    for (auto i = 0; i < l; i++) {
        g[i] = 1 << (torusBits - (i + 1) * radixBits); // 1/(B^(i) as Torus: 2^torusBits * 2^(-radixBits*i)
    }
    return g;
}

void gadgetDecompose(DecomposedData& out, const Integer in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    UnsignedInteger tmp = (out.sign == 1) ? in : -in;
    UnsignedInteger mask = ((1 << param.radixBits) - 1) << (param.torusBits - param.radixBits);
    for (auto i = 0; i < out.l; i++) {
        out.value[i] = (mask & tmp) >> (param.torusBits - (i + 1) * param.radixBits);
        mask >>= param.radixBits;
    }
}

void gadgetDecomposeNtt(DecomposedDataDft& out, const NttType in, const YatfheParameters& param) {
    NttType mask = ((1 << param.radixBits) - 1) << (param.dftBits - param.radixBits);
    for (auto i = 0; i < out.l; i++) {
        out.value[i] = (mask & in) >> (param.dftBits - (i + 1) * param.radixBits);
        mask >>= param.radixBits;
    }
}

Integer recompose(const DecomposedData& digits, const YatfheParameters& param) {
    Integer res {0};
    for (auto i = 0; i < digits.value.size(); ++i) {
        res += digits.value[i] << (param.torusBits - (i + 1) * param.radixBits);
    }
    return res * digits.sign;
}

/**
 * Calculate in * B^-j.
 * @param in The input to decompose. Bit length should be less than (maxIntegerBitLength - (torusBits - radixBits)).
 * i.e. (32 - (32 - 4)) = 4. Therefore, max in should be less than 2^3.
 * @param param
 * @return
 */
std::vector<Integer> decomposeOverB(const Integer in, const YatfheParameters& param) {
    std::vector<Integer> output(param.ksLevel);
    for (int i = 0; i < output.size(); ++i) {
        output[i] = in << (param.torusBits - (i + 1) * param.radixBits);
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
    copy(tmp.begin(), tmp.begin() + out.l, out.value.begin());
}

// todo: incorrect, need fix
void signedGadgetDecompositionNtt(DecomposedDataDft& out, const NttType in, const YatfheParameters& param) {
    vector<NttType> tmp(param.dftBits / param.radixBits);
    UnsignedInteger carry = 0;
    for (auto i = 0; i < tmp.size(); i++) {
        auto unsignedDigit = ((in >> (i * param.radixBits)) & param.digitMask) + carry;
        auto carryMask = unsignedDigit & param.baseOverTwo;
        auto signedDigit = unsignedDigit - (carryMask << 1);
        carry = carryMask >> (param.radixBits - 1);
        tmp[tmp.size() - i - 1] = signedDigit;
    }
    copy(tmp.begin(), tmp.begin() + out.l, out.value.begin());
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
void gadgetDecomposeTrlweNtt(DecomposedTrlwe& output, const TrlweDft& input, const YatfheParameters& param) {
    const auto k = input.k;
    const auto N = input.b.coeffs.size();
    const auto l = output.lDft;
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
void recomposeTrlweNtt(TrlweDft& output, const DecomposedTrlwe& input, const YatfheParameters& param) {
    const auto k = output.k;
    const auto N = output.b.coeffs.size();
    const auto l = input.lDft;
    trlweSetZero(output.a, output.b);
    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto row = 0; row < k + 1; row++) {
            auto& currIn = (row < k) ? input.rlweDfts[lvl].a[row] : input.rlweDfts[lvl].b;
            auto& currOut = (row < k) ? output.a[row] : output.b;
            for (auto j = 0; j < N; j++) {
                currOut.coeffs[j] = modAdd(currOut.coeffs[j], currIn.coeffs[j] << (64 - (lvl + 1) * param.radixBits));
            }
        }
    }
}

//// G^-1 * Trlwe = DecomposedTrlwe
//void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param) {
//    const auto k = input.k;
//    const auto N = input.b.coeffs.size();
//    const auto l = output.l;
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