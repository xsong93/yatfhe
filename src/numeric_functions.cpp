//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include "numeric_functions.h"
#include "torus.h"

using namespace std;
default_random_engine generator;
random_device rd;
mt19937 rng(rd());
uniform_int_distribution<Torus> uniformTorus32Distrib(INT32_MIN, INT32_MAX);

// Gaussian sample centered in message, with standard deviation sigma
Torus gaussian32(Torus message, const double sigma) {
    //Attention: all the implementation will use the stdev instead of the gaussian fourier param
    normal_distribution<double> distribution(0.,sigma); //TODO: can we create a global distrib of param 1 and multiply by sigma?
    double err = distribution(generator);
    return message + dToT32(err);
}

// Convert double to Torus32
Torus dToT32(const double d) {
    return int32_t(int64_t((d - int64_t(d)) * twoP32));
}