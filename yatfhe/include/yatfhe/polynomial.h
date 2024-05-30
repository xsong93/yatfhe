//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include "yatfhe/torus.h"

constexpr uint8_t POLY_MAX8 = 1 << 7;

struct TorusPolynomial {
    std::vector<Torus> coeffs {}; // N
    int N {};

    explicit TorusPolynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    TorusPolynomial(int N, Torus value) :
            coeffs(N, value),
            N(N) {};
};

using IntPolynomial = TorusPolynomial;
using BinPolynomial = IntPolynomial;

struct LagrangePolynomial {
    std::vector<NttType> coeffs {}; // N
    int N {};

    explicit LagrangePolynomial(int N) :
        N(N),
        coeffs(N, 0) {};

    LagrangePolynomial(int N, NttType value) :
        N(N),
        coeffs(N, value) {};
};

struct DoublePolynomial {
    std::vector<double> coeffs {}; // N
    int N {};

    explicit DoublePolynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    DoublePolynomial(int N, double value) :
            N(N),
            coeffs(N, value) {};
};

struct Ntt16Polynomial {
    std::vector<Ntt16> coeffs {}; // N
    int N {};

    explicit Ntt16Polynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    Ntt16Polynomial(int N, Ntt16 value) :
            N(N),
            coeffs(N, value) {};
};

struct Ntt14Polynomial {
    std::vector<Ntt14> coeffs {}; // N
    int N {};

    explicit Ntt14Polynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    Ntt14Polynomial(int N, Ntt16 value) :
            N(N),
            coeffs(N, value) {};
};

struct Ntt64Polynomial {
    std::vector<Ntt64> coeffs {};
    int N {};
    Ntt64Polynomial() : N(), coeffs() {};
    explicit Ntt64Polynomial (int n):
            coeffs(n,0), N(n) {};
};

void intPolyToDoublePoly(DoublePolynomial& output, const IntPolynomial & input);

void torusPolyToDoublePoly(DoublePolynomial& output, const TorusPolynomial& input);

void doublePolyToTorusPoly(TorusPolynomial& output, const DoublePolynomial& input);

void torusPolyToIntPoly(IntPolynomial& output, const TorusPolynomial& input, int mSize);

void intPolyToTorusPoly(TorusPolynomial& output, const IntPolynomial& input, int mSize);

void roundErrorPoly(DoublePolynomial& target, int torusBase);

void generateTestPolynomial(TorusPolynomial& v, int modP, int modQ);

void torusPolynomialRotate(TorusPolynomial& out, int a, const TorusPolynomial& input);

void torusPolynomialRotateMinusOne(TorusPolynomial& out, int a, const TorusPolynomial& input);

void polynomialMulNaive(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialMulAccNaive(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialAccumulate(TorusPolynomial& res, const TorusPolynomial& accum);

void polynomialAdd(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialAddSubOffset(TorusPolynomial& poly, int offset, bool isAdd);

void polynomialSub(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void generateLagrangePolynomialWithValueAt(LagrangePolynomial& lagrangePolynomial, int value, int position);

void lagrangePolynomialAccumulate(LagrangePolynomial& accum, LagrangePolynomial& poly);

void lagrangePolynomialAdd(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2);

void lagrangePolynomialSub(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2);

#endif //HLS_YATFHE_POLYNOMIAL_H
