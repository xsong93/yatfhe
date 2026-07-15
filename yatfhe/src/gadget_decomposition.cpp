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

// Signed (balanced) counterpart of gadgetDecomposeKs, over the KSK's own
// ksRadixBits/ksWidthBits. See signedGadgetDecomposition for the digit/residual
// rationale. switchKeyForTlwe consumes this as a signed MAC (coeff = value[j]*
// sign accumulated in int64, one longModP(., LWE_Q) at the end), so balanced
// digits drop in unchanged; the top carry-out has weight 2^ksWidthBits == LWE_Q
// and vanishes mod LWE_Q, keeping the key switch exact.
void signedGadgetDecompositionKs(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    const int radixBits = param.ksRadixBits;
    const int widthBits = param.ksWidthBits;
    const int l = out.l;
    const Torus B = static_cast<Torus>(1) << radixBits;
    const Torus halfB = B >> 1;

    out.sign = 1;

    // Low ksWidthBits bits (LWE_Q domain): two's complement of the signed input.
    UnsignedInteger u = static_cast<UnsignedInteger>(in);
    if (widthBits < static_cast<int>(sizeof(UnsignedInteger) * 8)) {
        u &= (static_cast<UnsignedInteger>(1) << widthBits) - 1;
    }

    // Round to the top l*ksRadixBits bits (residual in [-Delta/2, Delta/2)).
    const int shift = widthBits - l * radixBits;   // # bits dropped below the last kept digit
    if (shift > 0) {
        u += static_cast<UnsignedInteger>(1) << (shift - 1);
    }

    Torus carry = 0;
    for (int j = l - 1; j >= 0; --j) {
        const UnsignedInteger window =
            (u >> (widthBits - (j + 1) * radixBits)) & static_cast<UnsignedInteger>(B - 1);
        Torus digit = static_cast<Torus>(window) + carry;
        if (digit >= halfB) {          // >= B/2  ->  fold into [-B/2, B/2), carry up
            digit -= B;
            carry = 1;
        } else {
            carry = 0;
        }
        out.value[j] = digit;
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
 * Signed (balanced) gadget decomposition g^-1(x).
 *
 * Keeps the l most significant radix-B digits of x (B = 2^radixBits), aligned
 * with genGadgetVector's weights 2^(torusBits-(i+1)*radixBits) -- out.value[0]
 * is the most significant digit. Each digit is balanced into [-B/2, B/2) by
 * carrying into the next-more-significant digit, and the discarded low tail is
 * rounded to nearest so the gadget residual lies in [-Delta/2, Delta/2) with
 * Delta = 2^(torusBits-l*radixBits), instead of the [0, Delta) of the
 * sign-magnitude gadgetDecompose. This halves the digit spread (E[d^2]=B^2/12
 * vs B^2/3) and the residual (E[eps^2]=Delta^2/12 vs Delta^2/3).
 *
 * Digits carry their own sign, so out.sign stays +1 (recompose* multiply
 * value[i] by sign, so a global sign is unnecessary here). Negative torus
 * values decompose directly from their two's-complement bits (e.g. in = -1
 * yields all -1 digits), so no |in| split is needed.
 */
void signedGadgetDecomposition(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    const int radixBits = param.radixBits;
    const int torusBits = param.torusBits;
    const int l = out.l;
    const Torus B = static_cast<Torus>(1) << radixBits;
    const Torus halfB = B >> 1;

    out.sign = 1;

    UnsignedInteger u = static_cast<UnsignedInteger>(in);

    // Round to the top l*radixBits bits: add half a ULP at the truncation point
    // so the dropped tail becomes a nearest-rounding residual in [-Delta/2, Delta/2).
    const int shift = torusBits - l * radixBits;   // # bits dropped below the last kept digit
    if (shift > 0) {
        u += static_cast<UnsignedInteger>(1) << (shift - 1);
    }

    // Extract the l kept digits, least-significant first so the balancing carry
    // propagates toward the more-significant digit.
    Torus carry = 0;
    for (int j = l - 1; j >= 0; --j) {
        const UnsignedInteger window =
            (u >> (torusBits - (j + 1) * radixBits)) & static_cast<UnsignedInteger>(B - 1);
        Torus digit = static_cast<Torus>(window) + carry;
        if (digit >= halfB) {          // >= B/2  ->  fold into [-B/2, B/2), carry up
            digit -= B;
            carry = 1;
        } else {
            carry = 0;
        }
        out.value[j] = digit;
    }
}
