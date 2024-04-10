//
// Created by Xintong Song on 2024/4/8.
//

#include "yatfhe/ntt.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"

// todo
// https://www.nayuki.io/page/number-theoretic-transform-integer-dft

void applyNtt32(IntPolynomial& output, const IntPolynomial& input, const int root, const int mod) {
    const auto N = input.N;
    auto& out = output.coeffs;
    for (auto i = 0; i < N; i++) {
        Integer sum = 0;
        for (auto j = 0; j < N; j++) {
            auto k = (Integer)((long long)i * j % N);
            auto temp = (long long)input.coeffs[j] * modPow(root, k, mod) + sum;
            sum = (int)(temp % mod);
        }
        out[i] = sum;
    }
}

void applyIntt32(IntPolynomial& output, IntPolynomial& input, const int root, const int mod) {
    applyNtt32(output, input, reciprocal(root, mod), mod);
    auto scaler = reciprocal(input.N, mod);
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = (Integer)((long long)output.coeffs[i] * scaler % mod);
    }
}