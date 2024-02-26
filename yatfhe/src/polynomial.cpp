//
// Created by Xintong Song on 2023/12/25.
//
#include "iostream"
#include "polynomial.h"
#include "numeric_functions.h"

void initTorusPolynomial(TorusPolynomial& torusPolynomial, const int N) {
    torusPolynomial.coeffs.resize(N);
    torusPolynomial.N = N;
}

void initLagrangePolynomial(LagrangePolynomial& lagrangePolynomial, const int N) {
    lagrangePolynomial.coeffs.resize(N);
    lagrangePolynomial.N = N;
}

int validateRotator(const int a, const int N) {
    int aTrue = a % (2 * N);
    if (aTrue < 0) {
        aTrue += 2 * N;
    }
    return (aTrue < N) ? aTrue : aTrue - N;
}

// vj = ((pj / q) mod p) / p
void generateTestPolynomial(TorusPolynomial& v, const int modP, const int modQ) {
    const auto N = v.N;
    for (auto i = 0; i < N; i++) {
        int tmp = (modP * i / modQ) % modP;
        cout << (double) tmp / modP << " ";
        v.coeffs[i] = doubleToTorus32((double) tmp / modP);
    }
    cout << endl;
}

// output = (X^{a}) * input
// for example. N = 5
// if a = 3, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> x^3 * a0 + x^4 * a1 - a2 - x * a3 - x^2 * a4 -> (-a2, -a3, -a4, a0, a1)
// if a = -1, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> -x^4 * a0 + a1 + x * a2 + x^2 * a3 + x^3 * a4 -> (a1, a2, a3, a4, -a0) <=> a = 9 <=> -1 * (a = 4)
void torusPolynomialRotate(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    const auto aTrue = validateRotator(a, N);
    const auto isWrap = (aTrue < N) ? 1 : -1; // 1: no wrap; -1:wrap around
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
    }
}

// output = (X^{a} - 1) * input = x^a * input - input
void torusPolynomialRotateMinusOne(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    const auto aTrue = validateRotator(a, N);
    const auto isWrap = (aTrue < N) ? 1 : -1; // 1: no wrap; -1:wrap around
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