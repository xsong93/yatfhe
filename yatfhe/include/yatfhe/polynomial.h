//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include "yatfhe/torus.h"
#include "yatfhe/numeric_functions.h"

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

struct Int8Polynomial {
    std::vector<int8_t> coeffs {}; // N
    int N {};

    explicit Int8Polynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    Int8Polynomial(int N, int8_t value) :
            coeffs(N, value),
            N(N) {};
};

struct Int8PolynomialD {
    std::vector<std::vector<int8_t>> coeffs {}; // N * d
    int N {};
    int d {};

    Int8PolynomialD(int N, int d) :
            coeffs(N, std::vector<int8_t>(d, 0)),
            N(N),
            d(d){};
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

struct Ntt24Polynomial {
    std::vector<Ntt24> coeffs {}; // N
    int N {};

    explicit Ntt24Polynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    Ntt24Polynomial(int N, Ntt24 value) :
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

// res += accum
template<typename PolyType>
void polynomialAccumulateI32(PolyType& res, const PolyType& accum) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] += accum.coeffs[i];
    }
}

template<typename PolyType>
void polynomialAccumulateT32(PolyType& res, const PolyType& accum) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = modAddT32(res.coeffs[i], accum.coeffs[i]);
    }
}

template<typename PolyType, typename U>
void polynomialAccumulateModP(PolyType& res, const PolyType& accum, const U p) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = longModP(static_cast<int64_t>(res.coeffs[i]) + static_cast<int64_t>(accum.coeffs[i]), p);
    }
}

void intPolyToDoublePoly(DoublePolynomial& output, const IntPolynomial& input);

void torusPolyToDoublePoly(DoublePolynomial& output, const TorusPolynomial& input);

void doublePolyToTorusPoly(TorusPolynomial& output, const DoublePolynomial& input);

void torusPolyToIntPoly(IntPolynomial& output, const TorusPolynomial& input, int mSize);

void intPolyToTorusPoly(TorusPolynomial& output, const IntPolynomial& input, int mSize);

void roundErrorTorusPoly(TorusPolynomial& target, int torusBase);

void roundErrorDoublePoly(DoublePolynomial& target, int torusBase);

void generateTestPolynomial(TorusPolynomial& v, int modP, int modQ);

void torusPolynomialRotate(TorusPolynomial& out, int a, const TorusPolynomial& input);

void torusPolynomialRotateMinusOne(TorusPolynomial& out, int a, const TorusPolynomial& input);

void intPolynomialRotate(IntPolynomial& out, int a, const IntPolynomial& input, int64_t p);

void int8PolynomialRotate(Int8Polynomial& out, int a, const Int8Polynomial& input, int modP);

void int8PolynomialRotateMinusOne(Int8Polynomial& out, int a, const Int8Polynomial& input, int modP);

void polynomialMulNaiveT32(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialMulNaiveModQ(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2, int64_t q);

void polynomialMulNaiveI8(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, int q);

void polynomialMulAccNaiveI8(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, int q);

void polynomialMulAccNaiveI32(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

void polynomialMulAccNaiveT32(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialMulAccNaiveT32b(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

//void polynomialAccumulateI32(TorusPolynomial& res, const TorusPolynomial& accum);

void polynomialAddI32(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

void polynomialAddT32(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void polynomialAddSubOffset(IntPolynomial& poly, int offset, bool isAdd);

void polynomialSubI32(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

void polynomialSubT32(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void generateLagrangePolynomialWithValueAt(LagrangePolynomial& lagrangePolynomial, int value, int position);

void lagrangePolynomialAccumulate(LagrangePolynomial& accum, LagrangePolynomial& poly);

void lagrangePolynomialAdd(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2);

void lagrangePolynomialSub(LagrangePolynomial& output, const LagrangePolynomial& input1, const LagrangePolynomial& input2);

#endif //HLS_YATFHE_POLYNOMIAL_H
