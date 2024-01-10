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
    return message + dToT32(err);
}

// Convert double to Torus32
Torus dToT32(const double d) {
    return int32_t(int64_t((d - int64_t(d)) * twoP32));
}