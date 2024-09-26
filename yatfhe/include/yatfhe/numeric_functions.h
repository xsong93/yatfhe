//
// Created by Xintong Song on 2023/12/11.
//

#ifndef HLS_YATFHE_NUMERIC_FUNCTIONS_H
#define HLS_YATFHE_NUMERIC_FUNCTIONS_H

#include <cstdlib>
#include <random>
#include "yatfhe/torus.h"

using namespace std;

extern random_device rd;
extern mt19937 rng;
extern uniform_int_distribution<Binary> binaryDistrib;
extern uniform_int_distribution<Torus> uniformTorusDistrib;
static const int64_t twoP32 = INT64_C(1) << 32; // 2^32

int calLogBase2(int N);

Integer genIntUniformDist(int lowerBound, int upperBound);

uint64_t genUInt64UniformDist(uint64_t lowerBound, uint64_t upperBound);

Torus addGaussianNoise(Torus message, double sigma);

Torus modSwitchToTorus32(int32_t mu, uint32_t Msize);

int32_t modSwitchFromTorus32(Torus in, uint32_t newMod);

uint32_t modSwitchFromTorus32Pos(Torus in, uint32_t newMod);

Torus doubleToTorus32(double d);

double roundError(double in, int torusBase);

Torus roundTorusError(Torus in, int torusBase);

Integer roundErrorForShiftedTorus(Torus in, double sigma, int shift);

int intModP(int a, int p);

int64_t longModP(int64_t a, int64_t p);

int64_t barrettReduceT32(int64_t in);

Torus modAddT32(Torus in1, Torus in2);

Torus modSubT32(Torus in1, Torus in2);

Torus modMulT32(Torus in1, Torus in2);

int64_t modInverse(int64_t a, int64_t p);

double torus32ToDouble(Torus in);

void initCoeffsViaUniformDistribution(vector<Torus>& coeffs);

void initCoeffsWithGaussianNoiseSingleSample(vector<Torus>& coeffs, Torus msg, double sigma);

void initCoeffsWithGaussianNoiseMultiSample(std::vector<Torus>& coeffs, const std::vector<Torus>& msg, double sigma);

Integer modPow(Integer x, Integer y, Integer mod);

Integer reciprocal(Integer a, Integer mod);

bool isPrime(int num);

void uniquePrimeFactors(std::vector<int>& result, int n);

int sqrtFloor(int x);

template<typename T, typename R>
T modP(const T a, const R p) {
    auto b = a % p;
    if (b > p / 2 - 1) {
        b -= p;
    } else if (b < - p / 2) {
        b += p;
    }
    return b;
}

template<typename T, typename R>
T modMulQ(T in1, T in2, R q) {
    return static_cast<T>(modP(static_cast<R>(in1) * static_cast<R>(in2), q));
}

template <typename T>
void setCoeffsValue(vector<T> coeffs, T val) {
    for (auto i = 0; i < coeffs.size(); i++) {
        coeffs[i] = val;
    }
}

template <typename T, typename R>
void vectorMultConst(vector<T>& output, vector<T>& input1, R num) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = (T) (input1[i] * num);
    }
}

template <typename T, typename R>
void vectorDotMultConst(vector<T>& output, vector<T>& input1, vector<R>& nums) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = (T) (input1[i] * nums[i]);
    }
}

template <typename T>
void vectorAdd(vector<T>& output, vector<T>& input1, vector<T>& input2) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = input1[i] + input2[i];
    }
}

template <typename T>
void vectorSub(vector<T>& output, vector<T>& input1, vector<T>& input2) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = input1[i] - input2[i];
    }
}

#endif //HLS_YATFHE_NUMERIC_FUNCTIONS_H
