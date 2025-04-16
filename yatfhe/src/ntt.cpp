//
// Created by Xintong Song on 2024/1/11.
//
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/torus.h"

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
                temp_mult = fastmm_opt(output[j*block_size + k + gap],tw[i][tw_index]);
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
    output = input;
    Ntt64 temp_add, temp_sub, temp_mult;
//    Ntt64 tw_flag = 0;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (i + 1);
        block_size = 1 << (i + 1);
        gap = 1 << i;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
//                tw_flag = tw[i][tw_index];
                temp_add = modADDscale64(output[j * block_size + k], output[j * block_size + k + gap]);
                temp_sub = modSUBscale64(output[j * block_size + k], output[j * block_size + k + gap]);
                temp_mult = fastmm_opt(temp_sub, tw[i][tw_index]);
                output[j * block_size + k] = temp_add;
                output[j * block_size + k + gap] = temp_mult;
            }
        }
    }
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
    int64_t temp_ntt;
    int64_t temp_poly;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
            temp_poly = temp_ntt % TORUS_Q;
            if (temp_poly  < TORUS_MIN) {
                out.coeffs[i] = int32_t(temp_poly + TORUS_Q);
            } else if (temp_poly > TORUS_MAX){
                out.coeffs[i] = int32_t(temp_poly - TORUS_Q);
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

Ntt64 modMul(Ntt64 a, Ntt64 b) {
    mpz_t A, B, TEMP, P;
    mpz_init(A);
    mpz_init(B);
    mpz_init(TEMP);
    mpz_init(P);
    mpz_set_ui(A,a);
    mpz_set_ui(B,b);
    mpz_set_ui(P,MOD64);
    mpz_mul(TEMP, A, B);
    mpz_mod(TEMP, TEMP, P);
    return Ntt64(mpz_get_ui(TEMP));
}

// output_j = aj * bj mod p
void modularMult(std::vector<uint64_t>& output, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsB) {
    const auto N = output.size();
    for (auto j = 0; j < N; j++) {
        output[j] = fastmm_opt(coeffsA[j], coeffsB[j]);
    }
}

// b += a * s mod p
void modularAccumulate(std::vector<uint64_t>& res, const std::vector<uint64_t>& in1, const std::vector<uint64_t>& in2) {
    const auto N = res.size();
    Ntt64 tmp = 0;
    for (auto j = 0; j < N; j++) {
        if (in1[j] == 0 || in2[j] == 0) {
            continue;
        }
        tmp = fastmm_opt(in1[j], in2[j]);
        res[j] = modAdd(res[j], tmp);
    }
}

// res = aN * sN
void calModularInnerProductNtt(LagrangePolynomial& res, const vector<LagrangePolynomial>& in1, const vector<LagrangePolynomial>& in2) {
    for (auto i = 0; i < in1.size(); i++) {
        modularAccumulate(res.coeffs, in1[i].coeffs, in2[i].coeffs);
    }
}

void calModularInnerProductNtt(LagrangePolynomial& res, const LagrangePolynomial& in1, const LagrangePolynomial& in2) {
    modularAccumulate(res.coeffs, in1.coeffs, in2.coeffs);
}