//
// Created by Xintong Song on 2023/12/11.
//

#ifndef HLS_YATFHE_NUMERIC_FUNCTIONS_H
#define HLS_YATFHE_NUMERIC_FUNCTIONS_H

#include <cstdlib>
#include <random>
#include "yatfhe/torus.h"
#include "yatfhe/csprng.h"

using namespace std;

extern thread_local ChaCha20Rng rng;

extern uniform_int_distribution<Binary> binaryDistrib;

extern uniform_int_distribution<Integer> ternaryDistrib;

uniform_int_distribution<Torus> uniformTorusDistrib(Torus min, Torus max);

uniform_int_distribution<NttType> uniformNttDistrib(NttType min, NttType max);

int calLogBase2(int N);

Integer genIntUniformDist(Integer lowerBound, Integer upperBound);

uint64_t genUInt64UniformDist(uint64_t lowerBound, uint64_t upperBound);

Torus addTUniformNoise(Torus message, int b, const int64_t torusQ);

Torus modSwitchToTorusGeneral(int32_t mu, uint32_t mSize, int64_t torusQ);

int64_t modSwitchFromTorusGeneral(Torus in, int64_t newMod, int64_t torusQ);

Torus modSwitchToTorus32(int32_t mu, uint32_t mSize);

int32_t modSwitchFromTorus32(Torus in, uint32_t newMod);

Torus doubleToTorus32(double d);

double roundError(double in, int torusBase);

Torus roundTorusGeneralError(const Torus in, const int torusBase, const int64_t q);

Torus roundTorus32Error(const Torus in, const int torusBase);

Integer roundErrorForShiftedTorus(Torus in, double sigma, int shift);

int intModP(int a, int p);

int64_t longModP(int64_t a, int64_t p);

int64_t barrettReduceT32(int64_t in);

int64_t montgomoryReduceT32(int64_t in);

// Torus addTorus(Torus in1, Torus in2);

Torus subTorus(int64_t q, Torus in1, Torus in2);

Torus multTorus(int64_t q, Torus in1, Torus in2);

int64_t modInverse(int64_t a, int64_t p);

double torus32ToDouble(Torus in);

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, Torus min, Torus max);

void initNttCoeffsViaUniformDistribution(std::vector<NttType>& coeffs, NttType min, NttType max);

void initCoeffsWithTUniformNoiseSingleSample(vector<Torus>& coeffs, Torus msg, int pos, double sigma, const int64_t torusQ);

void initCoeffsWithTUniformNoiseMultiSample(std::vector<Torus>& coeffs, const std::vector<Torus>& msg, double sigma
                                            , const int64_t torusQ);

Integer modPow(Integer x, Integer y, Integer mod);

Integer reciprocal(Integer a, Integer mod);

bool isPrime(int num);

void uniquePrimeFactors(std::vector<int>& result, int n);

int sqrtFloor(int x);

template<typename T, typename R>
T modP(const T a, const R p) {
    auto b = a % p;
    if (b > (p - 1) / 2) {
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

template<typename... TorusArgs>
std::common_type_t<TorusArgs...> addTorus(int64_t q, TorusArgs... inputs) {
    if (q == Q_32) {
        return (inputs + ...);
    }
    int64_t tmp = (static_cast<int64_t>(inputs) + ...);
//    if (tmp > TORUS_MAX) {
//        return static_cast<std::common_type_t<TorusArgs...>>(tmp - TORUS_Q);
//    }
//    if (tmp < TORUS_MIN) {
//        return static_cast<std::common_type_t<TorusArgs...>>(tmp + TORUS_Q);
//    }
//    return static_cast<std::common_type_t<TorusArgs...>>(tmp);
    return longModP(tmp, q);
}

template <typename T, typename R>
void vectorMultConst(vector<T>& output, vector<T>& input1, R num) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = (T) (input1[i] * num);
    }
}

template <typename T, typename R, typename U>
void vectorMultConstModQ(vector<T>& output, vector<T>& input1, R num, U q) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = (T) modMulQ(input1[i], num, q);
    }
}

template <typename T, typename R>
void vectorDotMultConst(vector<T>& output, vector<T>& input1, vector<R>& nums) {
    for (auto i = 0; i < output.size(); i++) {
        output[i] = (T) (input1[i] * nums[i]);
    }
}

template <typename T>
void vectorAdd(vector<T>& output, vector<T>& input1, vector<T>& input2, int64_t q) {
    for (auto i = 0; i < output.size(); i++) {
//        output[i] = input1[i] + input2[i];
        output[i] = addTorus(q, input1[i], input2[i]);
    }
}

template <typename T>
void vectorSub(vector<T>& output, vector<T>& input1, vector<T>& input2, int64_t q) {
    for (auto i = 0; i < output.size(); i++) {
//        output[i] = input1[i] - input2[i];
        output[i] = subTorus(q, input1[i], input2[i]);
    }
}

#endif //HLS_YATFHE_NUMERIC_FUNCTIONS_H
