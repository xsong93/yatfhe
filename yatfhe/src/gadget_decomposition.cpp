//
// Created by Xintong Song on 2024/3/15.
//
#include <vector>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/gadget_decomposition.h"

using namespace std;

// offset = B/2 * (2^(torusBits - radixBits) + 2^(torusBits - 2 * radixBits) + ... + 2^(torusBits - l * radixBits))
int genOffset(const int radixBits, const int bHalf, const int l, const int torusBits) {
    int64_t res = 0;
    for (auto i = 1; i <= l; ++i) {
        res += static_cast<int64_t>(1) << (torusBits - i * radixBits);
    }
    return static_cast<int>(res * bHalf);
}

// g = (1/B, ..., 1/B^l), B = 2^radixBits
std::vector<Torus> genGadgetVector(const int radixBits, const int l, const int torusBits) {
    std::vector<Torus> g(l);
    for (auto i = 0; i < l; i++) {
        g[i] = static_cast<Torus>(1) << (torusBits - (i + 1) * radixBits); // 1/(B^(i) as Torus: 2^torusBits * 2^(-radixBits*i)
    }
    return g;
}

void gadgetDecompose(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    const int64_t in64 = static_cast<int64_t>(in);
    const uint64_t absIn = (in64 < 0) ? static_cast<uint64_t>(-in64) : static_cast<uint64_t>(in64);
    UnsignedInteger tmp = static_cast<UnsignedInteger>(absIn);
    UnsignedInteger mask = ((static_cast<UnsignedInteger>(1) << param.radixBits) - 1) << (param.torusBits - param.radixBits);
    for (auto i = 0; i < out.l; i++) {
        out.value[i] = static_cast<Torus>((mask & tmp) >> (param.torusBits - (i + 1) * param.radixBits));
        mask >>= param.radixBits;
    }
}

void gadgetDecomposeKs(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    const int64_t in64 = static_cast<int64_t>(in);
    const uint64_t absIn = (in64 < 0) ? static_cast<uint64_t>(-in64) : static_cast<uint64_t>(in64);
    UnsignedInteger tmp = static_cast<UnsignedInteger>(absIn);
    UnsignedInteger mask = ((static_cast<UnsignedInteger>(1) << param.ksRadixBits) - 1) << (param.ksWidthBits - param.ksRadixBits);
    for (auto i = 0; i < out.l; i++) {
        out.value[i] = static_cast<Torus>((mask & tmp) >> (param.ksWidthBits - (i + 1) * param.ksRadixBits));
        mask >>= param.ksRadixBits;
    }
}

Torus recomposeSelf(const DecomposedData& digits, const YatfheParameters& param) {
    Torus res {0};
    for (auto i = 0; i < digits.value.size(); ++i) {
        res += digits.value[i] << (param.torusBits - (i + 1) * param.radixBits);
    }
    return res * digits.sign;
}

void recomposeFirstHalf(DecomposedData& output, const DecomposedData& lhs, const std::vector<DecomposedData>& mid) {
    for (auto j = 0; j < lhs.l; j++) {
        for (auto l = 0; l < mid[0].l; l++) {
            output.value[j] += (lhs.value[l] * lhs.sign) * (mid[j].value[l] * mid[j].sign);
        }
    }
}

Torus recomposeTwoParts(const DecomposedData& lhs, const std::vector<Integer>& rhs) {
    Torus out {0};
    for (auto i = 0; i < rhs.size(); i++) {
        out += lhs.value[i] * rhs[i] * lhs.sign;
    }
    return out;
}

/**
 * Calculate in * B^-j.
 * @param in The input to decompose. Bit length should be less than (maxIntegerBitLength - (torusBits - radixBits)).
 * i.e. (32 - (32 - 4)) = 4. Therefore, max in should be less than 2^3.
 * @param param
 * @return
 */
void decomposeOverB(std::vector<Torus>& output, const Integer in, const YatfheParameters& param) {
    for (int i = 0; i < output.size(); ++i) {
        output[i] = static_cast<Torus>(in) << (param.ksWidthBits - (i + 1) * param.ksRadixBits);
    }
}

/**
 * Calculate signed decomposition g^-1(x).
 * @param res Resulting g^-1(x).
 * @param input The input to decompose.
 * @param param
 * //todo: need test
 */
void signedGadgetDecomposition(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    out.sign = (in < 0) ? -1 : 1;
    const int64_t in64 = static_cast<int64_t>(in);
    const uint64_t absIn = (in64 < 0) ? static_cast<uint64_t>(-in64) : static_cast<uint64_t>(in64);
    UnsignedInteger unsignedIn = static_cast<UnsignedInteger>(absIn);
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
