//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include <iostream>
#include "yautil/numeric_functions.h"
#include "yatfhe/torus.h"

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

Torus modSwitchToTorus32(int32_t mu, int32_t Msize) {
    uint64_t interv = ((UINT64_C(1) << 63) / Msize) * 2; // width of each interval
    uint64_t phase64 = mu * interv;
    //floor to the nearest multiples of interv
    return phase64 >> 32;
}

int32_t modSwitchFromTorus32(Torus phase, int32_t Msize) {
    uint64_t interv = ((UINT64_C(1) << 63) / Msize) * 2; // width of each interval
    uint64_t half_interval = interv / 2; // begin of the first intervall
    uint64_t phase64 = (uint64_t(phase) << 32) + half_interval;
    //floor to the nearest multiples of interv
    return (phase >= 0) ? (phase64 / interv) : (phase64 / interv - Msize);
}

// offset = B/2 * (2^(torusBits - radixBits) + 2^(torusBits - 2 * radixBits) + ... + 2^(torusBits - l * radixBits))
int genOffset(const int radixBits, const int bHalf, const int l, const int torusBits) {
    int res = 0;
    for (auto i = 1; i <= l; ++i) {
        res += 1 << (torusBits - i * radixBits);
    }
    return res * bHalf;
}

// g = (1/B, ..., 1/B^l), B = 2^radixBits
std::vector<Torus> genGadgetVector(const int radixBits, const int l, const int torusBits) {
    std::vector<Torus> g(l);
    for (auto i = 1; i <= l; i++) {
        g[i - 1] = 1 << (torusBits - i * radixBits); // 1/(B^(i) as Torus: 2^torusBits * 2^(-radixBits*i)
    }
    return g;
}

void signedGadgetDecomposition(vector<Torus>& res, const Torus input) {
    // todo
}

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, const int N) {
    for (int j = 0; j < N; j++) {
        coeffs[j] = uniformTorus32Distrib(rng);
    }
}

void initCoeffsWithGaussianNoise(std::vector<Torus>& coeffs, const Torus msg, const int N, const double sigma) {
    for (int j = 0; j < N; j++) {
        coeffs[j] = addGaussianNoise(msg, sigma);
    }
}