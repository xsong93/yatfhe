//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yatfhe/ntt.h"

int validateRotator(const int a, const int N) {
    int aTrue = a % (2 * N);
    if (aTrue < 0) {
        aTrue += 2 * N;
    }
    return (aTrue < N) ? aTrue : aTrue - N;
}

void intPolyToDoublePoly(DoublePolynomial& output, const IntPolynomial & input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = (double) input.coeffs[i];
    }
}

void torusPolyToDoublePoly(DoublePolynomial& output, const TorusPolynomial& input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = torus32ToDouble(input.coeffs[i]);
    }
}

void doublePolyToTorusPoly(TorusPolynomial& output, const DoublePolynomial& input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = doubleToTorus32(input.coeffs[i]);
    }
}

void torusPolyToIntPoly(IntPolynomial& output, const TorusPolynomial& input, const int mSize) {
    for (auto i = 0; i < input.N; i++) {
        output.coeffs[i] = modSwitchFromTorus32(input.coeffs[i], mSize);
    }
}

void intPolyToTorusPoly(TorusPolynomial& output, const IntPolynomial& input, const int mSize) {
    for (auto i = 0; i < input.N; i++) {
        output.coeffs[i] = modSwitchToTorus32(input.coeffs[i], mSize);
    }
}

void generateLagrangePolynomialWithValueAt(LagrangePolynomial& lagrangePolynomial, const int value, const int position) {
    IntPolynomial tmp(lagrangePolynomial.N);
    tmp.coeffs[position] = value;
    applyNtt(lagrangePolynomial, tmp);
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

// accum += poly
void lagrangePolynomialAccumulate(LagrangePolynomial& accum, LagrangePolynomial& poly) {
    const auto N = accum.N;
    for (auto i = 0; i < N; i++) {
        accum.coeffs[i] += modAdd(accum.coeffs[i], poly.coeffs[i]);
    }
}

// res = poly1 + poly2
void polynomialAdd(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] + poly2.coeffs[i];
    }
}

/**
 * Add or sub a value to every coefficients of the target polynomial.
 * @param poly Target polynomial.
 * @param offset Offset value.
 * @param isAdd True: add offset. Otherwise, subtract offset.
 */
void polynomialAddSubOffset(TorusPolynomial& poly, const int offset, const bool isAdd) {
    for (auto i = 0; i < poly.N; i++) {
        poly.coeffs[i] = isAdd ? (poly.coeffs[i] + offset) : (poly.coeffs[i] - offset);
    }
}

// res = poly1 - poly2
void polynomialSub(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] - poly2.coeffs[i];
    }
}
