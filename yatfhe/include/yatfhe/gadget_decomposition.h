//
// Created by Xintong Song on 2024/3/15.
//

#ifndef YATFHE_GADGET_DECOMPOSITION_H
#define YATFHE_GADGET_DECOMPOSITION_H

#include <limits>
#include <stdexcept>
#include <string>
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

// a window is the radixBits bits below its shift
constexpr inline int windowShift(const int widthBits, const int radixBits, const int level) {
    return widthBits - (level + 1) * radixBits;
}

// Rounding bit of the tail below the last kept digit: adding it rounds that tail half away from zero.
// Zero when the l*radixBits kept bits are all the bits there are.
constexpr inline UnsignedInteger roundingBias(const int radixBits, const int widthBits, const int l) {
    const int droppedBits = widthBits - l * radixBits;
    return droppedBits > 0 ? static_cast<UnsignedInteger>(1) << (droppedBits - 1) : 0;
}

// Truncates a value to widthBits bits.
// A torus narrower than UnsignedInteger keeps its two's complement sign.
constexpr inline UnsignedInteger maskToWidth(const UnsignedInteger value, const int widthBits) {
    constexpr int widthOfUnsigned = std::numeric_limits<UnsignedInteger>::digits;
    return widthBits >= widthOfUnsigned ? value : value & ((static_cast<UnsignedInteger>(1) << widthBits) - 1);
}

constexpr inline Torus balancedDigit(const Torus digit, const Torus radixBase, const bool carries) {
    return carries ? digit - radixBase : digit;
}

// Rounds then decomposes one value, handing emit(level, digit) from level l-1 (least significant) down to level 0.
template <typename Emit>
inline void zeroMeanDigits(UnsignedInteger u, const int l, const int radixBits, const int widthBits, Emit&& emit) {
    const UnsignedInteger widthMask = maskToWidth(~UnsignedInteger{0}, widthBits);
    u &= widthMask;
    const bool negative = (u >> (widthBits - 1)) & 1;
    const UnsignedInteger m = (negative ? (UnsignedInteger{0} - u) & widthMask : u) + roundingBias(radixBits, widthBits, l);

    const Torus kRadixBase = static_cast<Torus>(1) << radixBits;
    const Torus kHalfB = kRadixBase >> 1;
    const Torus sign = negative ? -1 : 1;

    Torus carry = 0;
    for (int level = l - 1; level >= 0; --level) {
        const int shift = windowShift(widthBits, radixBits, level);
        const auto window = static_cast<Torus>((m >> shift) & static_cast<UnsignedInteger>(kRadixBase - 1));
        const Torus digit = window + carry;
        const Torus tieBit = level > 0 ? static_cast<Torus>((m >> (shift + 2 * radixBits - 1)) & 1) : 0;
        carry = (digit > kHalfB) | ((digit == kHalfB) & tieBit);
        emit(level, sign * balancedDigit(digit, kRadixBase, carry));
    }
}

inline void zeroMeanDigits(DecomposedData& out, const Torus in, const int radixBits, const int widthBits) {
    zeroMeanDigits(static_cast<UnsignedInteger>(in), out.l, radixBits, widthBits, [&out](const int level, const Torus digit){
        out.value[level] = digit;
    });
}

inline void signedGadgetDecomposition(DecomposedData& out, const Torus in, const YatfheParameters& param) {
    out.sign = 1;
    zeroMeanDigits(out, in, param.radixBits, param.torusBits);
}

// Decomposes a whole row, writing digit j of every value to level j, one level at a time.
// The per-value state is held in flat arrays, which is what lets the inner loop vectorise.
template <typename OutT>
void decomposeRowByLevel(OutT* const* outPtr, const Torus* in, const int N, const int L, const int radixBits,
                         const int widthBits) {
    const Torus kRadixBase = static_cast<Torus>(1) << radixBits;
    const Torus kHalfB = kRadixBase >> 1;
    const UnsignedInteger widthMask = maskToWidth(~UnsignedInteger{0}, widthBits);
    const UnsignedInteger round = roundingBias(radixBits, widthBits, L);

    thread_local std::vector<UnsignedInteger> magnitudeBuf;
    thread_local std::vector<Torus> signBuf;
    thread_local std::vector<Torus> carryBuf;
    if (magnitudeBuf.size() < static_cast<size_t>(N)) {
        magnitudeBuf.resize(N);
        signBuf.resize(N);
        carryBuf.resize(N);
    }
    UnsignedInteger* __restrict magnitude = magnitudeBuf.data();
    Torus* __restrict sign = signBuf.data();
    Torus* __restrict carry = carryBuf.data();

    for (int j = 0; j < N; j++) {
        const UnsignedInteger value = static_cast<UnsignedInteger>(in[j]) & widthMask;
        const bool negative = (value >> (widthBits - 1)) & 1; // top bit: the sign in two's complement
        magnitude[j] = (negative ? (UnsignedInteger{0} - value) & widthMask : value) + round;
        sign[j] = negative ? -1 : 1;
        carry[j] = 0;
    }

    for (int level = L - 1; level >= 0; --level) {
        const int shift = windowShift(widthBits, radixBits, level);
        // no window above the top digit, so the tie read is masked off.
        const int tieShift = level > 0 ? shift + 2 * radixBits - 1 : 0;
        const Torus tieMask = level > 0 ? 1 : 0;
        OutT* __restrict out = outPtr[level];
        for (int j = 0; j < N; j++) {
            const auto window = static_cast<Torus>((magnitude[j] >> shift) & static_cast<UnsignedInteger>(kRadixBase - 1));
            const Torus digit = window + carry[j];
            const Torus tieBit = static_cast<Torus>((magnitude[j] >> tieShift) & 1) & tieMask;
            carry[j] = (digit > kHalfB) | ((digit == kHalfB) & tieBit);
            out[j] = static_cast<OutT>((digit - carry[j] * kRadixBase) * sign[j]);
        }
    }
}

// Row-wise gadget decomposition, same digits as signedGadgetDecomposition
template <int L, typename OutT>
void decomposeRowUnrolled(OutT* const* outPtr, const Torus* in, const int N, const int radixBits, const int widthBits) {
    const Torus kRadixBase = static_cast<Torus>(1) << radixBits;
    const Torus kHalfB = kRadixBase >> 1;
    const UnsignedInteger widthMask = maskToWidth(~UnsignedInteger{0}, widthBits);
    const UnsignedInteger round = roundingBias(radixBits, widthBits, L);

    if constexpr (L < 5) {
        // Coefficient-major: one pass per value keeps its whole digit sequence in registers.
        OutT* levelOut[L];
        for (int level = 0; level < L; level++) {
            levelOut[level] = outPtr[level];
        }
        for (int j = 0; j < N; j++) {
            const UnsignedInteger value = static_cast<UnsignedInteger>(in[j]) & widthMask;
            const bool negative = (value >> (widthBits - 1)) & 1; // top bit: the sign in two's complement
            const UnsignedInteger magnitude = (negative ? (UnsignedInteger{0} - value) & widthMask : value) + round;

            Torus carry = 0;
            for (int level = L - 1; level >= 0; --level) {
                const int shift = windowShift(widthBits, radixBits, level);
                const auto window = static_cast<Torus>((magnitude >> shift) & static_cast<UnsignedInteger>(kRadixBase - 1));
                const Torus digit = window + carry;
                const Torus tieBit = level > 0 ? static_cast<Torus>((magnitude >> (shift + 2 * radixBits - 1)) & 1) : 0;
                carry = (digit > kHalfB) | ((digit == kHalfB) & tieBit);
                const Torus balanced = balancedDigit(digit, kRadixBase, carry);
                levelOut[level][j] = static_cast<OutT>(negative ? -balanced : balanced);
            }
        }
    } else {
        decomposeRowByLevel(outPtr, in, N, L, radixBits, widthBits);
    }
}

template<typename OutT>
void decomposeRow(OutT* const* outPtr, const Torus* in, const int N, const int l, const YatfheParameters& param) {
    const int radixBits = param.radixBits;
    const int widthBits = param.torusBits;

    if (l * radixBits > widthBits) {
        throw std::invalid_argument("decomposeRow: " + std::to_string(l) +
                                        " digits of " + std::to_string(radixBits) +
                                        " bits do not fit the " + std::to_string(widthBits) + "-bit torus");
    }
    if constexpr (std::numeric_limits<OutT>::is_signed) {
        // digits reach +2^(radixBits-1), which needs radixBits - 1 < the type's value bits
        if (radixBits - 1 >= std::numeric_limits<OutT>::digits) {
            throw std::invalid_argument("decomposeRow: digits up to +2^" + std::to_string(radixBits - 1) +
                                        " do not fit the " + std::to_string(std::numeric_limits<OutT>::digits + 1) + "-bit digit type."
                                        " Use radixBits <= " + std::to_string(std::numeric_limits<OutT>::digits));
        }
    }

    switch (l) {
        case 1:
            decomposeRowUnrolled<1>(outPtr, in, N, radixBits, widthBits);
            return;
        case 2:
            decomposeRowUnrolled<2>(outPtr, in, N, radixBits, widthBits);
            return;
        case 3:
            decomposeRowUnrolled<3>(outPtr, in, N, radixBits, widthBits);
            return;
        case 4:
            decomposeRowUnrolled<4>(outPtr, in, N, radixBits, widthBits);
            return;
        case 5:
            decomposeRowUnrolled<5>(outPtr, in, N, radixBits, widthBits);
            return;
        case 6:
            decomposeRowUnrolled<6>(outPtr, in, N, radixBits, widthBits);
            return;
        case 7:
            decomposeRowUnrolled<7>(outPtr, in, N, radixBits, widthBits);
            return;
        case 8:
            decomposeRowUnrolled<8>(outPtr, in, N, radixBits, widthBits);
            return;
        default:
            decomposeRowByLevel(outPtr, in, N, l, radixBits, widthBits);
    }
}

int genOffset(int radixBits, int bHalf, int l, int torusBits);

#endif //YATFHE_GADGET_DECOMPOSITION_H
