//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include "torus.h"

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

void initTorusPolynomial(TorusPolynomial& torusPolynomial, int N);

void initLagrangePolynomial(LagrangePolynomial& lagrangePolynomial, int N);

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

void deletePolynomial(TorusPolynomial& polynomial);

void deletePolynomial(LagrangePolynomial& polynomial);

#endif //HLS_YATFHE_POLYNOMIAL_H
