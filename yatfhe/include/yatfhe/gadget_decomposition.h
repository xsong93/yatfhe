//
// Created by Xintong Song on 2024/3/15.
//

#ifndef YATFHE_GADGET_DECOMPOSITION_H
#define YATFHE_GADGET_DECOMPOSITION_H

#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"

struct DecomposedData {
    std::vector<Torus> value; // l
    int l {};
    Integer sign {1}; // set default to 1 as positive sign

    explicit DecomposedData(int size) : l(size), value(size, 0) {};
};

struct DecomposedDataDft {
    std::vector<NttType> value; // l
    int l {};

    explicit DecomposedDataDft(int size) : l(size), value(size, 0) {};
};

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

void gadgetDecompose(DecomposedData& out, Torus in, const YatfheParameters& param);

// Same as gadgetDecompose but uses the KSK's own ksRadixBits/ksWidthBits
// instead of PBS's radixBits/torusBits. Only for key switch.
void gadgetDecomposeKs(DecomposedData& out, Torus in, const YatfheParameters& param);

// Only for key switch.
void signedGadgetDecompositionKs(DecomposedData& out, Torus in, const YatfheParameters& param);

Torus recomposeSelf(const DecomposedData& digits, const YatfheParameters& param);

void recomposeFirstHalf(DecomposedData& output, const DecomposedData& lhs, const std::vector<DecomposedData>& mid);

Torus recomposeTwoParts(const DecomposedData& lhs, const std::vector<Integer>& rhs);

void decomposeOverB(std::vector<Torus>& output, Integer in, int bitWidth, int radixBits);

void decomposeOverBKS(std::vector<Torus>& output, Integer in, const YatfheParameters& param);

/**
 * Signed gadget decomposition.
 * Each digit is balanced into [-B/2, B/2) by carrying into the next-more-significant digit, and the discarded low tail is
 * rounded to nearest so the gadget residual lies in [-Delta/2, Delta/2), with Delta = 2^(torusBits-l*radixBits).
 */
inline void signedGadgetDecomposition(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    const int radixBits = param.radixBits;
    const int torusBits = param.torusBits;
    const int l = out.l;
    const Torus B = static_cast<Torus>(1) << radixBits;
    const Torus halfB = B >> 1;

    out.sign = 1;

    auto u = static_cast<UnsignedInteger>(in);

    // Round to the top l*radixBits bits
    const int shift = torusBits - l * radixBits; // # bits dropped below the last kept digit
    if (shift > 0) {
        u += static_cast<UnsignedInteger>(1) << (shift - 1);
    }

    // Extract the l kept digits, least-significant first so the balancing carry propagates toward the more-significant digit.
    Torus carry = 0;
    for (int j = l - 1; j >= 0; --j) {
        const UnsignedInteger window = (u >> (torusBits - (j + 1) * radixBits)) & static_cast<UnsignedInteger>(B - 1);
        Torus digit = static_cast<Torus>(window) + carry;
        carry = digit >= halfB; // >= B/2  ->  fold into [-B/2, B/2), carry up
        digit -= carry * B;
        out.value[j] = digit;
    }
}

// Row-wise gadget decomposition
template <int L, typename OutT>
void decomposeRowUnrolled(OutT* const* outPtr, const Torus* in, const int N, const int radixBits, const int torusBits) {
    // Coefficient-major
    // Threshold at 5 keeps the fast path for both.
    const Torus B = static_cast<Torus>(1) << radixBits;
    const Torus halfB = B >> 1;
    const int shift = torusBits - L * radixBits;
    const UnsignedInteger round = shift > 0 ? static_cast<UnsignedInteger>(1) << (shift - 1) : 0;

    if constexpr (L < 5) {
        OutT* o[L];
        for (int lvl = 0; lvl < L; lvl++) {
            o[lvl] = outPtr[lvl];
        }
        for (int j = 0; j < N; j++) {
            const UnsignedInteger u = static_cast<UnsignedInteger>(in[j]) + round;
            Torus carry = 0;
            for (int lvl = L - 1; lvl >= 0; --lvl) {
                const UnsignedInteger window = (u >> (torusBits - (lvl + 1) * radixBits)) & static_cast<UnsignedInteger>(B - 1);
                Torus digit = static_cast<Torus>(window) + carry;
                carry = digit >= halfB;
                digit -= carry * B;
                o[lvl][j] = static_cast<OutT>(digit);
            }
        }
    } else {
        // digit-major
        thread_local std::vector<UnsignedInteger> uBuf;
        thread_local std::vector<Torus> carryBuf;
        if (uBuf.size() < static_cast<size_t>(N)) {
            uBuf.resize(N);
            carryBuf.resize(N);
        }
        UnsignedInteger* __restrict u = uBuf.data();
        Torus* __restrict c = carryBuf.data();
        for (int j = 0; j < N; j++) {
            u[j] = static_cast<UnsignedInteger>(in[j]) + round;
            c[j] = 0;
        }
        for (int lvl = L - 1; lvl >= 0; --lvl) {
            const int sh = torusBits - (lvl + 1) * radixBits;
            OutT* __restrict out = outPtr[lvl];
            for (int j = 0; j < N; j++) {
                const auto w = static_cast<Torus>((u[j] >> sh) & static_cast<UnsignedInteger>(B - 1));
                const Torus digit = w + c[j];
                const Torus carry = digit >= halfB ? 1 : 0;
                c[j] = carry;
                out[j] = static_cast<OutT>(digit - carry * B);
            }
        }
    }
}

template<typename OutT>
void decomposeRow(OutT* const* outPtr, const Torus* in, const int N, const int l, const YatfheParameters& param) {
    const int b = param.radixBits, t = param.torusBits;
    switch (l) {
        case 1: decomposeRowUnrolled<1>(outPtr, in, N, b, t); return;
        case 2: decomposeRowUnrolled<2>(outPtr, in, N, b, t); return;
        case 3: decomposeRowUnrolled<3>(outPtr, in, N, b, t); return;
        case 4: decomposeRowUnrolled<4>(outPtr, in, N, b, t); return;
        case 5: decomposeRowUnrolled<5>(outPtr, in, N, b, t); return;
        case 6: decomposeRowUnrolled<6>(outPtr, in, N, b, t); return;
        case 7: decomposeRowUnrolled<7>(outPtr, in, N, b, t); return;
        case 8: decomposeRowUnrolled<8>(outPtr, in, N, b, t); return;
        default: {
            DecomposedData d{l};
            for (auto j = 0; j < N; j++) {
                signedGadgetDecomposition(d, in[j], param);
                for (auto lvl = 0; lvl < l; lvl++) {
                    outPtr[lvl][j] = static_cast<OutT>(d.value[lvl] * d.sign);
                }
            }
        }
    }
}

int genOffset(int radixBits, int bHalf, int l, int torusBits);

#endif //YATFHE_GADGET_DECOMPOSITION_H
