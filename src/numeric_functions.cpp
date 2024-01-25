//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include "numeric_functions.h"
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