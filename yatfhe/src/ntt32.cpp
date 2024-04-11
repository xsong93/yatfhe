//
// Created by Xintong Song on 2024/4/8.
//

#include "yatfhe/ntt.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"

// todo: support negative input, support negacyclic property
// https://www.nayuki.io/page/number-theoretic-transform-integer-dft

int findModulus(int vecLen, int minimum) {
    int start = std::max((minimum - 1 + vecLen - 1) / vecLen, 1);
    for (long n = (long)start * vecLen + 1; n <= std::numeric_limits<int>::max(); n += vecLen) {
        if (isPrime((int)n)) {
            return (int)n;
        }
    }
    throw std::runtime_error("Modulus not found");
}

void applyNtt32(IntPolynomial& output, const IntPolynomial& input, const int root, const int mod) {
    const auto N = input.N;
    auto& out = output.coeffs;
    for (auto i = 0; i < N; i++) {
        int sum = 0;
        for (auto j = 0; j < N; j++) {
            auto k = (int)((long long)i * j % N);
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
        output.coeffs[i] = (int)((long long)output.coeffs[i] * scaler % mod);
    }
}

int modMul32(const int in1,const int in2, const int mod) {
    return (in1 * in2) % mod;
}

void circularConvolve(IntPolynomial& output, IntPolynomial& poly1, IntPolynomial& poly2, const int mod) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = modMul32(poly1.coeffs[i], poly2.coeffs[i], mod);
    }
}

bool isPrimitiveRoot(int val, int degree, int mod) {
    if (val < 0 || val >= mod) {
        throw std::invalid_argument("Invalid value");
    }
    if (degree < 1 || degree >= mod) {
        throw std::invalid_argument("Invalid degree");
    }

    if (modPow(val, degree, mod) != 1) {
        return false;
    }

    std::vector<int> primeFactors;
    uniquePrimeFactors(primeFactors, degree);
    for (int p : primeFactors) {
        if (modPow(val, degree / p, mod) == 1)
            return false;
    }
    return true;
}

int findPrimitiveRoot(int degree, int totient, int mod) {
    if (degree < 1 || degree > totient || totient >= mod || totient % degree != 0) {
        throw std::invalid_argument("Invalid arguments");
    }

    for (int g = 1; g < mod; g++) {
        if (isPrimitiveRoot(g, totient, mod)) {
            return modPow(g, totient / degree, mod);
        }
    }

    throw std::runtime_error("Primitive root not found");
}