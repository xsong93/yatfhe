//
// Created by Xintong Song on 2023/12/11.
//

#ifndef HLS_YATFHE_NUMERIC_FUNCTIONS_H
#define HLS_YATFHE_NUMERIC_FUNCTIONS_H

#include <cstdlib>
#include <random>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/yatfhe_parameters.h"

using namespace std;

extern random_device rd;
extern mt19937 rng;
extern uniform_int_distribution<Torus> uniformTorus32Distrib;
static const int64_t twoP32 = INT64_C(1) << 32; // 2^32

Torus addGaussianNoise(Torus message, double sigma);

int32_t modSwitchFromTorus32(Torus phase, int32_t Msize);

Torus doubleToTorus32(double d);

double torus32ToDouble(Torus in);

Torus modSwitchToTorus32(int32_t mu, int32_t Msize);

Torus int2torus(uint64_t x, int log_scale);

int32_t genOffset(int bgBit, int halfBg, int l);

void gadgetDecomposition(vector<vector<IntPolynomial>>& output, const vector<TorusPolynomial>& input, const YatfheParameters& param);

void modularAccumulate(vector<uint64_t>& coeffsB, const vector<uint64_t>& coeffsA, const vector<uint64_t>& coeffsS);

void calModularInnerProductNtt(LagrangePolynomial& b, LagrangePolynomial& a, const IntPolynomial& s, const int N);

void initCoeffsViaUniformDistribution(vector<Torus>& coeffs, int N);

void initCoeffsWithGaussianNoise(vector<Torus>& coeffs, Torus msg, int N, double sigma);

#endif //HLS_YATFHE_NUMERIC_FUNCTIONS_H
