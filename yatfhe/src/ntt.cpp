//
// Created by Xintong Song on 2024/1/11.
//
#include <iostream>
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"

using namespace std;

void DITNRLaPoly(LagrangePolynomial& out, const LagrangePolynomial& in) {
    auto& output = out.coeffs;
    const auto& input = in.coeffs;
    const auto& N = in.N;
    output = input;

    const auto& tw = NWC_TW64.tw_factor;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt64 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
        for (auto j = 0; j < block; j++) {
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                temp_mult = modMULT(output[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modAdd(output[j*block_size + k], temp_mult);
                temp_sub = modSub(output[j*block_size + k], temp_mult);

                output[j*block_size + k] = temp_add;
                output[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void DIFRNLaPoly(LagrangePolynomial& out, const LagrangePolynomial& in) {
    auto &output = out.coeffs;
    const auto &input = in.coeffs;
    const auto &tw = NWC_ITW64.tw_factor;
    const auto &N = in.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
//    Ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    output = input;
    Ntt64 temp_add, temp_sub, temp_mult;
//    int pos_a, pos_b;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (i + 1);
        block_size = 1 << (i + 1);
        gap = 1 << i;
//        pos_a = 0;
//        pos_b = 0;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
//                flag_a = res[j * block_size + k];
//                flag_b = res[j * block_size + k + gap];
//                flag_tw = tw[i][tw_index];
//                pos_a = j * block_size + k;
//                pos_b = j * block_size + k + gap;
                temp_add = modADDscale(output[j * block_size + k], output[j * block_size + k + gap]);
                temp_sub = modSUBscale(output[j * block_size + k], output[j * block_size + k + gap]);
                temp_mult = modMULT(temp_sub, tw[i][tw_index]);
                output[j * block_size + k] = temp_add;
                output[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void initGlobalParamsNtt64(int N) {
    auto depth = clog2(N);
    TwParam::initTwParam(NWC_TW64, depth);
    TwParam::initTwParam(NWC_ITW64, depth);
    TwRom::initTwRom(TW_ROM64, N);
    genTW_ROM(TW_ROM64);
    genNWCparam(NWC_TW64, N, TW_ROM64, STR_NTT);
    genNWCparam(NWC_ITW64, N, TW_ROM64, STR_INTT);
}

// Function to perform Number Theoretic Transform (NTT)
void applyNtt(LagrangePolynomial& out, const IntPolynomial& in) {
    auto N = in.N;
    LagrangePolynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (in.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt64(in.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt64(in.coeffs[i] + MOD64);
        }
    }
    DITNRLaPoly(out, format_input);
}

void applyIntt(IntPolynomial& out, const LagrangePolynomial& in) {
    auto N = in.N;
    LagrangePolynomial res(N);
    DIFRNLaPoly(res, in);
    int64_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT64_MASK);
        if (temp_poly >= POLY_MAX) {
            out.coeffs[i] = int32_t(temp_poly - (POLY_MAX<<1));
        } else {
            out.coeffs[i] = int32_t(temp_poly);
        }
    }
}

void bitRevShuffle(std::vector<NttType>& x) {
    int j = 0;
    int b = 0;
    int N = int(x.size());
    for (int i = 1; i < N; i++) {
        b = N >> 1;  // Initialize b to half of N
        while (j >= b) {
            j -= b;  // Perform bit-reversal
            b >>= 1;
        }
        j += b;  // Move to the next position

        // Swap elements if the bit-reversed index is greater than the current index
        if (j > i) {
            NttType temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

Ntt64 modAdd(Ntt64 x, Ntt64 y) {
    return ((MOD64 - x) > y) ? (x + y) : (x + y - MOD64);
}

Ntt64 modSub(Ntt64 x, Ntt64 y) {
    return (x >= y) ? (x - y) : (MOD64 - y + x);
}

// todo: debug
Ntt64 modMul(Ntt64 x, Ntt64 y) {
    // Break down x and y into 32-bit components
    auto x0 = (uint32_t)x;
    auto x1 = (uint32_t)(x >> 32);
    auto y0 = (uint32_t)y;
    auto y1 = (uint32_t)(y >> 32);

    // Perform 64-bit multiplication
    uint64_t x0y0 = (uint64_t)x0 * (uint64_t)y0;
    uint64_t x0y1 = (uint64_t)x0 * (uint64_t)y1;
    uint64_t x1y0 = (uint64_t)x1 * (uint64_t)y0;
    uint64_t x1y1 = (uint64_t)x1 * (uint64_t)y1;

    // Compute partial products and handle carry
    auto d = (uint32_t)x0y0;
    uint64_t pp1 = (x0y0 >> 32) + (uint32_t)(x1y0) + (uint32_t)(x0y1);
    auto c = (uint32_t)pp1;
    uint64_t pp2 = (x1y0 >> 32) + (x0y1 >> 32) + (uint32_t)(x1y1);
    uint64_t pp3 = (pp1 >> 32) + (uint32_t)(pp2);

    // Handle overflow and underflow
    uint32_t a = (pp2 >> 32) + (x1y1 >> 32);
    uint64_t bpc = (uint32_t)pp3 + (uint64_t)c;
    bpc = ((bpc + (bpc >> 32)) << 32) - (bpc >> 32);
    uint64_t minus = ((uint64_t)a + ((uint64_t)(uint32_t)pp3));
    uint64_t plus = bpc + (uint64_t)d;

    // Return the result modulo MODULUS
    if (plus >= minus) {
        return (plus - minus);
    }
    return MOD64 - minus + plus;
}

// output_j = aj * bj mod p
void modularMult(std::vector<uint64_t>& output, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsB) {
    const auto N = output.size();
    for (auto j = 0; j < N; j++) {
        output[j] = modMULT(coeffsA[j], coeffsB[j]);
    }
}

// b += a * s mod p
void modularAccumulate(std::vector<uint64_t>& coeffsB, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsS) {
    const auto N = coeffsB.size();
    for (auto j = 0; j < N; j++) {
        auto tmp = modMULT(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modAdd(coeffsB[j], tmp);
    }
}

// b = aN * sN
void calModularInnerProductNtt(LagrangePolynomial& b, const vector<LagrangePolynomial>& a, const vector<LagrangePolynomial>& s) {
    for (auto i = 0; i < a.size(); i++) {
        modularAccumulate(b.coeffs, a[i].coeffs, s[i].coeffs);
    }
}