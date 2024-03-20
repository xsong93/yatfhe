//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include <iostream>
#include <iomanip>
#include "yatfhe/numeric_functions.h"
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"

using namespace std;
random_device rd;
mt19937 rng(rd());
uniform_int_distribution<Torus> uniformTorus32Distrib(INT32_MIN, INT32_MAX);

// Gaussian sample centered in message, with standard deviation sigma
Torus addGaussianNoise(Torus message, const double sigma) {
    normal_distribution<double> normalDistribution(0.0, sigma);
    double err = normalDistribution(rng);
    return message + doubleToTorus32(err);
}

// Convert double to Torus32, d in [-0.5, 0.5)
Torus doubleToTorus32(const double d) {
    auto frac = d - (int64_t) d; // get the fraction part of d
    if (frac >= 0.5) {
        frac = frac - 1;
    } else if (frac < -0.5) {
        frac = 1 + frac;
    }
    auto scaledFrac = int64_t(frac * twoP32); // scale the fraction part to [0, 2^32), then cast the result to a 64-bit integer
    return int32_t(scaledFrac); // rescale to int32
}

double torus32ToDouble(const Torus in) {
    return double(in) / twoP32;
}

double roundError(const double in, const int torusBase) {
    int mulP  = round(in * torusBase);
    int modP = mulP % torusBase;
    printf("yatfhe/src/numeric_functions.cpp@roundError. in: %f, mulP: %d, modP: %d\n", in, mulP, modP);
    return modP / double(torusBase);
}

Torus modSwitchToTorus32(int32_t mu, int32_t Msize) {
    uint64_t interv = ((UINT64_C(1) << 63) / Msize) * 2; // width of each interval
    uint64_t phase64 = mu * interv;
    //floor to the nearest multiples of interv
    return phase64 >> 32;
}

int32_t modSwitchFromTorus32(Torus in, int32_t newMod) {
    uint64_t interv = ((UINT64_C(1) << 63) / newMod) * 2; // width of each interval
    uint64_t half_interval = interv / 2; // begin of the first intervall
    uint64_t phase64 = (uint64_t(in) << 32) + half_interval;
    //floor to the nearest multiples of interv
    return (in >= 0) ? (phase64 / interv) : (phase64 / interv - newMod);
}

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs) {
    for (int& coeff : coeffs) {
        coeff = uniformTorus32Distrib(rng);
    }
}

void initCoeffsWithGaussianNoise(std::vector<Torus>& coeffs, const Torus msg, const double sigma) {
    for (int& coeff : coeffs) {
        coeff = addGaussianNoise(msg, sigma);
    }
}