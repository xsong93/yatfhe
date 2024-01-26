//
// Created by Xintong Song on 2023/12/25.
//
#include "polynomial.h"

void initTorusPolynomial(TorusPolynomial& torusPolynomial, const int N) {
    torusPolynomial.coeffs.resize(N);
    torusPolynomial.N = N;
}

void initLagrangePolynomial(LagrangePolynomial& lagrangePolynomial, const int N) {
    lagrangePolynomial.coeffs.resize(N);
    lagrangePolynomial.N = N;
}

// output = (X^{a} - 1) * input = x^a * input - input
void torusPolynomialMulByXaiMinusOne(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const int32_t N = input.N;
    const int32_t aTrue = (a < N) ? a : a - N;
    for (int32_t i = 0; i < N; i++) {
        out.coeffs[i] = ((i < aTrue) ? -input.coeffs[i - aTrue + N] : input.coeffs[i - aTrue]) - input.coeffs[i];
    }
}

void deletePolynomial(TorusPolynomial& polynomial) {
//    delete polynomial.coeffs;
//    polynomial.coeffs = nullptr;
}

void deletePolynomial(LagrangePolynomial& polynomial) {
//    delete polynomial.coeffs;
//    polynomial.coeffs = nullptr;
}