//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_H
#define HLS_YATFHE_NTT_H

#include <vector>
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"

constexpr uint32_t POLY_MAX = 1 << 31;
constexpr uint64_t POLY_Q = 1l << 32;
const string STR_NTT = "NWC-DIT-NR-NNT";
const string STR_INTT = "NWC-DIF-RN-INNT";
using namespace std;

int findModulus(int vecLen, int minimum);

bool isPrimitiveRoot(int g, int degree, int mod);

int findPrimitiveRoot(int degree, int totient, int mod);

void applyNtt(LagrangePolynomial& out, const IntPolynomial& in);

void applyNtt32(IntPolynomial& output, const IntPolynomial& input, int root, int mod);

void applyIntt32(IntPolynomial& output, IntPolynomial& input, int root, int mod);

int modMul32(int in1, int in2, int mod);

void circularConvolve(IntPolynomial& output, IntPolynomial& poly1, IntPolynomial& poly2, int mod);

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

template <typename T>
void bitRev(std::vector<T>& x) {
    int j = 0;
    int b = 0;
    int N = int(x.size());
    for (int i = 1; i < N; i++) {
        b = N >> 1;  // Initialize b to half of N
        while (j >= b) {
            j -= b;  // Perform bit-reversal
            b >>= 1;
        }
        j += b;  // Move to the next position

        // Swap elements if the bit-reversed index is greater than the current index
        if (j > i) {
            T temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

void bitRevShuffle(std::vector<NttType>& x);

NttType modAdd(NttType x, NttType y);

NttType modSub(NttType x, NttType y);

NttType modMul(NttType x, NttType y);

void modularMult(std::vector<NttType>& output, const std::vector<NttType>& coeffsA, const std::vector<NttType>& coeffsB);

void modularAccumulate(vector<NttType>& res, const vector<NttType>& in1, const vector<NttType>& in2);

void calModularInnerProductNtt(LagrangePolynomial& res, const vector<LagrangePolynomial>& in1, const vector<LagrangePolynomial>& in2);

void calModularInnerProductNtt(LagrangePolynomial& res, const LagrangePolynomial& in1, const LagrangePolynomial& in2);

#endif //HLS_YATFHE_NTT_H
