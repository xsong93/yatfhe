//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_OLD_H
#define HLS_YATFHE_NTT_OLD_H

#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"

constexpr NttType MODULUS = 0xffffffff00000001UL;
using namespace std;

void applyNttOld(NttPolynomial& out, const IntPolynomial& in);

void applyNtt32Old(IntPolynomial& output, const IntPolynomial& input, int root, int mod);

void applyIntt32Old(IntPolynomial& output, IntPolynomial& input, int root, int mod);

void applyNttTorusOld(NttPolynomial& out, const TorusPolynomial & in, int mSize);

void applyInttOld(IntPolynomial& out, const NttPolynomial& in);

void bitRevShuffleOld(std::vector<NttType>& x, int N);

NttType modAddOld(NttType x, NttType y);

NttType modSubOld(NttType x, NttType y);

NttType modMulOld(NttType x, NttType y);

void modularMultOld(std::vector<NttType>& output, const std::vector<NttType>& coeffsA, const std::vector<NttType>& coeffsB);

void modularAccumulateOld(vector<NttType>& coeffsB, const vector<NttType>& coeffsA, const vector<NttType>& coeffsS);

void calModularInnerProductNttOld(NttPolynomial& b, const vector<NttPolynomial>& a, const vector<NttPolynomial>& s);

#endif //HLS_YATFHE_NTT_OLD_H
