//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include "yatfhe/torus.h"

struct LagrangePolynomial {
    std::vector<uint64_t> coeffs {}; // N
    int N {};

    explicit LagrangePolynomial(int N) :
        N(N),
        coeffs(N, 0) {};

    explicit LagrangePolynomial(int N, uint64_t value) :
        N(N),
        coeffs(N, value) {};
};

struct DoublePolynomial {
    std::vector<double> coeffs {}; // N
    int N {};

    explicit DoublePolynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    explicit DoublePolynomial(int N, uint64_t value) :
            N(N),
            coeffs(N, value) {};
};

void intPolyToDoublePoly(DoublePolynomial& output, const IntPolynomial & input);

void torusPolyToDoublePoly(DoublePolynomial& output, const TorusPolynomial& input);

void doublePolyToTorusPoly(TorusPolynomial& output, const DoublePolynomial& input);

void torusPolyToIntPoly(IntPolynomial& output, const TorusPolynomial& input, int mSize);

void intPolyToTorusPoly(TorusPolynomial& output, const IntPolynomial& input, int mSize);

void generateLagrangePolynomialWithValueAt(LagrangePolynomial& lagrangePolynomial, int value, int position);

void generateTestPolynomial(TorusPolynomial& v, int modP, int modQ);

void torusPolynomialRotate(TorusPolynomial& out, int a, const TorusPolynomial& input);

void torusPolynomialRotateMinusOne(TorusPolynomial& out, int a, const TorusPolynomial& input);

void polynomialMulNaive(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialAccumulate(TorusPolynomial& res, const TorusPolynomial& accum);

void lagrangePolynomialAccumulate(LagrangePolynomial& accum, LagrangePolynomial& poly);

void polynomialAdd(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialAddSubOffset(TorusPolynomial& poly, int offset, bool isAdd);

void polynomialSub(TorusPolynomial& res, const IntPolynomial& poly1, const TorusPolynomial& poly2);

#endif //HLS_YATFHE_POLYNOMIAL_H
