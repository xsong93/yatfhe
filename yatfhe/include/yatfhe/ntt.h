//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_H
#define HLS_YATFHE_NTT_H

#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"

constexpr uint64_t MODULUS = 0xffffffff00000001UL;
using namespace std;

void applyNtt(LagrangePolynomial& out, const IntPolynomial& in);

void applyNttTorus(LagrangePolynomial& out, const TorusPolynomial & in, int mSize);

void applyIntt(IntPolynomial& out, LagrangePolynomial& in);

template <typename T, typename R>
void applyNttForAB(T& out, R& in) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyNtt(out.a[row], in.a[row]);
    }
    applyNtt(out.b, in.b);
}

template <typename T, typename R>
void applyInttForAB(T& out, R& in) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyIntt(out.a[row], in.a[row]);
    }
    applyIntt(out.b, in.b);
}

void bitRevShuffle(std::vector<uint64_t>& x, int N);

uint64_t modAdd(uint64_t x, uint64_t y);

uint64_t modSub(uint64_t x, uint64_t y);

uint64_t modMul(uint64_t x, uint64_t y);

void modularMult(std::vector<uint64_t>& output, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsB);

void modularAccumulate(vector<uint64_t>& coeffsB, const vector<uint64_t>& coeffsA, const vector<uint64_t>& coeffsS);

void calModularInnerProductNtt(LagrangePolynomial& b, LagrangePolynomial& a, const LagrangePolynomial& s);

#endif //HLS_YATFHE_NTT_H
