//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"

constexpr uint8_t POLY_MAX8 = 1 << 7;

struct TorusPolynomial {
    std::vector<Torus> coeffs {}; // N
    int N {};

    TorusPolynomial() = default;

    explicit TorusPolynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    TorusPolynomial(int N, Torus value) :
            coeffs(N, value),
            N(N) {};

    int minIndex() const {
        const auto it = std::min_element(coeffs.begin(), coeffs.end());
        return static_cast<int>(std::distance(coeffs.begin(), it));
    }

    int maxIndex() const {
        const auto it = std::max_element(coeffs.begin(), coeffs.end());
        return static_cast<int>(std::distance(coeffs.begin(), it));
    }
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

struct Int16Polynomial {
    std::vector<int16_t> coeffs {}; // N
    int N {};

    explicit Int16Polynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    Int16Polynomial(int N, int16_t value) :
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

struct DecompPolynomial {
    std::vector<Decomp> coeffs {}; // N
    int N {};

    DecompPolynomial() = default;

    explicit DecompPolynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    DecompPolynomial(int N, Decomp value) :
            coeffs(N, value),
            N(N) {};
};

using IntPolynomial = TorusPolynomial;
using BinPolynomial = IntPolynomial;

struct NttPolynomial {
    std::vector<NttType> coeffs {}; // N
    int N {};

    NttPolynomial() = default;

    explicit NttPolynomial(int N) :
        N(N),
        coeffs(N, 0) {};

    NttPolynomial(int N, NttType value) :
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

struct Ntt32Polynomial {
    std::vector<Ntt32> coeffs {}; // N
    int N {};

    explicit Ntt32Polynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    Ntt32Polynomial(int N, Ntt32 value) :
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

template<typename... PolyArgs>
void addTorusPolynomial(TorusPolynomial& res, const PolyArgs&... polys) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        int64_t tmp = (polys.coeffs[i] + ...);
        res.coeffs[i] = longModP(tmp, TORUS_Q);
    }
}

template<typename PolyType>
void accumulateTorusPolynomial(PolyType& res, const PolyType& accum) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = addTorus(TORUS_Q, res.coeffs[i], accum.coeffs[i]);
    }
}

template<typename PolyType, typename U>
void accumulatePolynomialModP(PolyType& res, const PolyType& accum, const U p) {
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

void generateTestPolynomialLt1(TorusPolynomial& v, int t);

void generateTestPolynomialCompLeq0(TorusPolynomial& v, const int t);

void generateTestPolynomialCompWithValue(TorusPolynomial& tv, int t, Integer v);

void generateTestPolynomialOne(TorusPolynomial& v);

void generateTestPolynomialValue(TorusPolynomial& tv, Integer v);

void validateRotator(int& aTrue, int& isWrap, int a, int N);

void rotateTorusPolynomial(TorusPolynomial& out, int a, const TorusPolynomial& input);

void rotateAccumulateTorusPolynomial(TorusPolynomial& accum, int aTrue, int isWrap, const TorusPolynomial& input);

void rotateTorusPolynomialMinusOne(TorusPolynomial& out, int a, const TorusPolynomial& input);

void rotateIntPolynomial(IntPolynomial& out, int a, const IntPolynomial& input, int64_t p);

void rotateInt8Polynomial(Int8Polynomial& out, int a, const Int8Polynomial& input, int modP);

void rotateInt8PolynomialMinusOne(Int8Polynomial& out, int a, const Int8Polynomial& input, int modP);

void multTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void multIntPolynomialModQ(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2, int64_t q);

void multInt8Polynomial(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, int q);

void multInt8PolynomialAcc(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, int q);

void multIntPolynomialAcc(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

void multTorusPolynomialAcc(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void addIntPolynomial(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

// void addTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void addSubIntPolynomialWithOffset(IntPolynomial& poly, int offset, bool isAdd);

void subIntPolynomial(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2);

void subTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2);

void rotateNttPolynomial(NttPolynomial& res, const NttPolynomial& in, int r);

void rotateNttPolynomialMinusOne(NttPolynomial& res, const NttPolynomial& in, int r);

void genNttPolynomialWithValueAt(NttPolynomial& lagrangePolynomial, int value, int position);

void accumulateNttPolynomial(NttPolynomial& accum, NttPolynomial& poly);

// void addNttPolynomial(NttPolynomial& output, const NttPolynomial& input1, const NttPolynomial& input2);

void subNttPolynomial(NttPolynomial& output, const NttPolynomial& input1, const NttPolynomial& input2);

void inverseGadgetDecomposePolynomial(vector<TorusPolynomial>& output, const IntPolynomial& input, const YatfheParameters& param);

#endif //HLS_YATFHE_POLYNOMIAL_H
