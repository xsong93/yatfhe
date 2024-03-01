//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include <iostream>
#include "numeric_functions.h"
#include "ntt.h"
#include "torus.h"

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

// offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
int32_t genOffset(const int bgBit, const int halfBg, const int l) {
    int32_t temp1 = 0;
    for (int32_t i = 0; i < l; ++i) {
        int32_t temp0 = 1 << (32 - (i + 1) * bgBit);
        temp1 += temp0;
    }
    return temp1 * halfBg;
}

// 1/B, ..., 1/B^l, B = 2^b
std::vector<Torus> genPowersOfBgbit(const int bgBit, const int l) {
    std::vector<Torus> h(l);
    for (auto i = 0; i < l; i++) {
        int power = (32 - (i + 1) * bgBit);
        h[i] = 1 << power; // 1/(bg^(i + 1)) as Torus32: 2^32 * 2^(-b*(i+1))
    }
    return h;
}

// output_j = aj * bj mod p
void modularMult(std::vector<uint64_t>& output, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsB) {
    const auto N = output.size();
    for (auto j = 0; j < N; j++) {
        output[j] = modMul(coeffsA[j], coeffsB[j]);
    }
}

// b += a * s mod p
void modularAccumulate(std::vector<uint64_t>& coeffsB, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsS) {
    const auto N = coeffsB.size();
    for (auto j = 0; j < N; j++) {
        auto tmp = modMul(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modAdd(coeffsB[j], tmp);
    }
}

// b = aN * sN
void calModularInnerProductNtt(LagrangePolynomial& b, LagrangePolynomial& a, const LagrangePolynomial& s) {
//    LagrangePolynomial sDft {N};
//    applyNtt(sDft, s);
    modularAccumulate(b.coeffs, a.coeffs, s.coeffs);
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