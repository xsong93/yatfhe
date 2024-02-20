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

// output = (X^{a}) * input
void torusPolynomialRotate(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    const auto aTrue = (a < N) ? a : a - N;
    const auto isWrap = (a < N) ? 1 : -1; // 1: no wrap; -1:wrap around
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
    }
}

// output = (X^{a} - 1) * input = x^a * input - input
void torusPolynomialRotateMinusOne(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    const auto aTrue = (a < N) ? a : a - N;
    const auto isWrap = (a < N) ? 1 : -1; // 1: no wrap; -1:wrap around
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = ((i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap)) - input.coeffs[i];
    }
}

void polynomialMulNaive(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        Torus tmp {0};
        for (int j = 0; j < N; j++) {
            tmp = (j <= i) ? (tmp + poly1.coeffs[j] * poly2.coeffs[i - j]) : (tmp - poly1.coeffs[j] * poly2.coeffs[N + i - j]);
        }
        res.coeffs[i] = tmp;
    }
}

// res += accum
void polynomialAccumulate(TorusPolynomial& res, const TorusPolynomial& accum) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] += accum.coeffs[i];
    }
}

// res = poly1 + poly2
void polynomialAdd(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] + poly2.coeffs[i];
    }
}

// res = poly1 - poly2
void polynomialSub(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] - poly2.coeffs[i];
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