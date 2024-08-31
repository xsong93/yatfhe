//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yautil/tool.h"

/**
 * For a random rotator input, this method converts the rotator to a value within the range of polynomial length.
 * Besides, in order to correctly show the negacyclic property, use a bit indicator to keep track of the negative signs after rotation.
 * For example. N = 5
 * if a = 3, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> x^3 * a0 + x^4 * a1 - a2 - x * a3 - x^2 * a4 -> (-a2, -a3, -a4, a0, a1)
 * if a = -1, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> -x^4 * a0 + a1 + x * a2 + x^2 * a3 + x^3 * a4 -> (a1, a2, a3, a4, -a0) <=> a = 9 <=> -1 * (a = 4)
 * @param aTrue Minimized rotator a.
 * @param isWrap 1: no wrap. -1: wrap around.
 * @param a Original rotator.
 * @param N Polynomial length,
 */
void validateRotator(int& aTrue, int& isWrap, const int a, const int N) {
    aTrue = a % (2 * N);
    if (aTrue < 0) {
        aTrue += 2 * N;
    }
    isWrap = (aTrue < N) ? 1 : -1;
    aTrue = (aTrue < N) ? aTrue : aTrue - N;
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

void roundErrorPoly(DoublePolynomial& target, const int torusBase) {
    for (auto i = 0 ; i < target.N; i++) {
        target.coeffs[i] = roundError(target.coeffs[i], torusBase);
    }
}

// vj = ((pj / q) mod p) / p
void generateTestPolynomial(TorusPolynomial& v, const int modP, const int modQ) {
    for (auto i = 0; i < v.N; i++) {
        int tmp = intModP((int)std::round(modP * i / modQ), modP);
        v.coeffs[i] = modSwitchToTorus32(tmp, modP);
    }
}

// output = (X^{a}) * input
void torusPolynomialRotate(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
    }
}

// output = (X^{a} - 1) * input = x^a * input - input
void torusPolynomialRotateMinusOne(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = ((i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap)) - input.coeffs[i];
    }
}

void int8PolynomialRotate(Int8Polynomial& out, const int a, const Int8Polynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * (int8_t)isWrap) : (input.coeffs[i - aTrue] * (int8_t)isWrap);
    }
}

void int8PolynomialRotateMinusOne(Int8Polynomial& out, const int a, const Int8Polynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    for (auto i = 0; i < N; i++) {
        out.coeffs[i] = ((i < aTrue) ? (-input.coeffs[i - aTrue + N] * (int8_t)isWrap) : (input.coeffs[i - aTrue] * (int8_t)isWrap)) - input.coeffs[i];
    }
}

void polynomialMulNaive(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    Torus tmp;
    for (auto i = 0; i < N; i++) {
        tmp = 0;
        for (auto j = 0; j < N; j++) {
            tmp = (j <= i) ? (tmp + poly1.coeffs[j] * poly2.coeffs[i - j]) : (tmp - poly1.coeffs[j] * poly2.coeffs[N + i - j]);
        }
        res.coeffs[i] = tmp;
    }
}

void polynomialMulAccNaive(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    Torus tmp;
    for (auto i = 0; i < N; i++) {
        tmp = 0;
        for (auto j = 0; j < N; j++) {
            tmp = (j <= i) ? (tmp + poly1.coeffs[j] * poly2.coeffs[i - j]) : (tmp - poly1.coeffs[j] * poly2.coeffs[N + i - j]);
        }
        res.coeffs[i] += tmp;
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
void polynomialAdd(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
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
void polynomialSub(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] - poly2.coeffs[i];
    }
}

void generateLagrangePolynomialWithValueAt(LagrangePolynomial& lagrangePolynomial, const int value, const int position) {
    IntPolynomial tmp {lagrangePolynomial.N};
    tmp.coeffs[position] = value;
    applyNtt(lagrangePolynomial, tmp);
}

// accum += poly
void lagrangePolynomialAccumulate(LagrangePolynomial& accum, LagrangePolynomial& poly) {
    const auto N = accum.N;
    for (auto i = 0; i < N; i++) {
        accum.coeffs[i] += modAdd(accum.coeffs[i], poly.coeffs[i]);
    }
}

void lagrangePolynomialAdd(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2) {
    for (auto i = 0; i < input1.N; i++) {
        output.coeffs[i] = modAdd(input1.coeffs[i], input2.coeffs[i]);
    }
}

void lagrangePolynomialSub(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2) {
    for (auto i = 0; i < input1.N; i++) {
        output.coeffs[i] = modSub(input1.coeffs[i], input2.coeffs[i]);
    }
}
