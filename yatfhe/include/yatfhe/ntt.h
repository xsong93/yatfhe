//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_H
#define HLS_YATFHE_NTT_H

#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"

constexpr NttType MODULUS = 0xffffffff00000001UL;
using namespace std;

void applyNtt(LagrangePolynomial& out, const IntPolynomial& in);

void applyNtt32(IntPolynomial& output, const IntPolynomial& input, int root, int mod);

void applyIntt32(IntPolynomial& output, IntPolynomial& input, int root, int mod);

void applyNttTorus(LagrangePolynomial& out, const TorusPolynomial & in, int mSize);

void applyIntt(IntPolynomial& out, const LagrangePolynomial& in);

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

void bitRevShuffle(std::vector<NttType>& x, int N);

NttType modAdd(NttType x, NttType y);

NttType modSub(NttType x, NttType y);

NttType modMul(NttType x, NttType y);

void modularMult(std::vector<NttType>& output, const std::vector<NttType>& coeffsA, const std::vector<NttType>& coeffsB);

void modularAccumulate(vector<NttType>& coeffsB, const vector<NttType>& coeffsA, const vector<NttType>& coeffsS);

void calModularInnerProductNtt(LagrangePolynomial& b, const vector<LagrangePolynomial>& a, const vector<LagrangePolynomial>& s);

#endif //HLS_YATFHE_NTT_H
