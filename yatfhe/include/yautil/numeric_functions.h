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

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

int genOffset(const int radixBits, const int bHalf, const int l, const int torusBits);

void initCoeffsViaUniformDistribution(vector<Torus>& coeffs, int N);

void initCoeffsWithGaussianNoise(vector<Torus>& coeffs, Torus msg, int N, double sigma);

#endif //HLS_YATFHE_NUMERIC_FUNCTIONS_H
