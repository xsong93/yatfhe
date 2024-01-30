//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
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

// Convert double to Torus32
Torus doubleToTorus32(const double d) {
    return int32_t(int64_t((d - int64_t(d)) * twoP32));
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
    return phase64 / interv;
}

//void modSwitchFromTorus32(int32_t& res, Torus phase, int32_t Msize) {
//    uint64_t interv = ((UINT64_C(1) << 63) / Msize) * 2; // width of each interval
//    uint64_t half_interval = interv / 2; // begin of the first intervall
//    uint64_t phase64 = (uint64_t(phase) << 32) + half_interval;
//    //floor to the nearest multiples of interv
//    res = phase64 / interv;
//}

Torus int2torus(uint64_t x, int log_scale) {
    const uint64_t bit_size = sizeof(Torus) * 8;
    return x << (bit_size - log_scale);
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

void gadgetDecomposition(vector<vector<IntPolynomial>>& output, const vector<TorusPolynomial>& input, const YatfheParameters& param) {
    const int k = param.k;
    const int N = param.N;
    const int l = param.l;
    const int bgBit = param.bgBit;
    const int maskMod = param.maskMod;
    const int halfBg = param.halfBg;
    const int offset = genOffset(bgBit, halfBg, l);
    vector<TorusPolynomial> buffer(input.size(), TorusPolynomial(param.N, 0));
    for (int i = 0; i <= k; i++) {
        for (int j = 0; j < N; j++) {
            buffer[i].coeffs[j] = input[i].coeffs[j] + offset;
        }
        for (int p = 0; p < l; p++) {
            const int decal = (32 - (p + 1) * bgBit);
            for (int j = 0; j < N; j++) {
                int32_t temp = (buffer[i].coeffs[j] >> decal) & maskMod;
                output[i][p].coeffs[j] = temp - halfBg;
            }
        }
    }
}

// b += a * s mod p
void modularAccumulate(std::vector<uint64_t>& coeffsB, std::vector<uint64_t>& coeffsA, std::vector<uint64_t>& coeffsS, const int N) {
    for (int j = 0; j < N; j++) {
        auto tmp = modMul(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modAdd(coeffsB[j], tmp);
    }
}

// b = akN * skN
void calModularInnerProduct(LagrangePolynomial& b, std::vector<TorusPolynomial>& a, const std::vector<IntPolynomial>& s, const int N, const int k) {
    for (int i = 0; i < k; i++) {
        initCoeffsViaUniformDistribution(a[i].coeffs, N);
        LagrangePolynomial sDft {N};
        LagrangePolynomial aDft {N};
        applyNtt(sDft, s[i]);
        applyNtt(aDft, a[i]);
        modularAccumulate(b.coeffs, aDft.coeffs, sDft.coeffs, N);
    }
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