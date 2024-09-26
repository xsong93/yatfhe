//
// Created by Xintong Song on 2023/12/25.
//
#include <random>
#include <iostream>
#include "yatfhe/numeric_functions.h"
#include "yatfhe/torus.h"

using namespace std;
random_device rd;
mt19937 rng(rd());
uniform_int_distribution<Binary> binaryDistrib(0, 1);
uniform_int_distribution<Torus> uniformTorusDistrib(TORUS_MIN, TORUS_MAX);

int calLogBase2(int N) {
    int res = 0;
    while (N >>= 1){
        res ++;
    }
    return res;
}

Integer genIntUniformDist(const int lowerBound, const int upperBound) {
    uniform_int_distribution<Integer> uniformIntDistrib(lowerBound, upperBound);
    return uniformIntDistrib(rng);
}

uint64_t genUInt64UniformDist(const uint64_t lowerBound, const uint64_t upperBound) {
    uniform_int_distribution<uint64_t> uniformIntDistrib(lowerBound, upperBound);
    return uniformIntDistrib(rng);
}

// Gaussian sample centered in message, with standard deviation sigma
Torus addGaussianNoise(Torus message, const double sigma) {
    normal_distribution<double> normalDistribution(0.0, sigma);
    Torus err = doubleToTorus32(normalDistribution(rng));
    Torus tmp = message + err;
    if ((message > 0 && tmp < 0) || (message < 0 && tmp > 0)) { // handle overflow
        return message - err;
    }
    return message + err;
}

/*
 Convert double to Torus32, d in [-0.5, 0.5).
 Generally, d shouldn't be taken too close to the boundary values, since in lwe-like cryptology schemes, a random error
 will be added for security considerations. An input value near boundary will very likely produce the add zero overflow issue.
*/
Torus doubleToTorus32(const double d) {
    auto frac = d - static_cast<double>(static_cast<int64_t>(d)); // get the fraction part of d
    if (frac >= 0.5) {
        frac = frac - 1;
    } else if (frac < -0.5) {
        frac = 1 + frac;
    }
    return static_cast<Torus>(frac * static_cast<double>(TORUS_Q)); // scale the fraction part to Z_TORUS_Q
}

double torus32ToDouble(const Torus in) {
    return static_cast<double>(in) / static_cast<double>(TORUS_Q);
}

double roundError(const double in, const int torusBase) {
    int mulP  = round(in * torusBase);
    int modP = mulP % torusBase;
//    printf("yatfhe/src/numeric_functions.cpp@roundError. in: %f, mulP: %d, modP: %d\n", in, mulP, modP);
    return modP / double(torusBase);
}

Torus roundTorusError(const Torus in, const int torusBase) {
    auto t = modSwitchFromTorus32(in, torusBase);
    auto t2 = modSwitchToTorus32(t, torusBase);
    return t2;
}

// To compensate the possible negative gaussian error, add sufficient sigma to the value so that it stops at the
// positive side of the nearest desired shift scale. The original integer value can then be restored after a proper right
// shift operation.
Integer roundErrorForShiftedTorus(const Torus in, const double sigma, const int shift) {
    return (in + 20 * doubleToTorus32(sigma)) >> shift;
}

int intModP(const int a, const int p) {
    auto b = a % p;
    if (b > p / 2 - 1) {
        b -= p;
    } else if (b < - p / 2) {
        b += p;
    }
    return b;
}

int64_t longModP(const int64_t a, const int64_t p) {
    auto b = a % p;
    if (b > p / 2 - 1) {
        b -= p;
    } else if (b < - p / 2) {
        b += p;
    }
    return b;
}

Torus modAddT32(Torus in1, Torus in2) {
    auto tmp = static_cast<int64_t>(in1) + static_cast<int64_t>(in2);
    return (tmp > TORUS_MAX) ? static_cast<Torus>(tmp - TORUS_Q) : static_cast<Torus>((tmp < TORUS_MIN) ? (TORUS_Q + tmp) : tmp);
}

Torus modSubT32(Torus in1, Torus in2) {
    auto tmp = static_cast<int64_t>(in1) - static_cast<int64_t>(in2);
    return (tmp > TORUS_MAX) ? static_cast<Torus>(tmp - TORUS_Q) : static_cast<Torus>((tmp < TORUS_MIN) ? (TORUS_Q + tmp) : tmp);
}

// Multiplicative inverse modulo p
int64_t modInverse(int64_t a, int64_t mod) {
    long m0 = mod, t, q;
    long x0 = 0, x1 = 1;

    if (mod == 1) {
        return 0;
    }

    while (a > 1) {
        q = a / mod; // q is quotient
        t = mod;
        mod = a % mod; // m is remainder now, process same as Euclid's algorithm
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

//    // Make x1 positive
//    if (x1 < 0) {
//        x1 += m0;
//    }

    return x1;
}

Torus modSwitchToTorus32(int32_t mu, uint32_t Msize) {
//    uint64_t interv = ((UINT64_C(1) << 63) / Msize) * 2; // width of each interval
//    uint64_t phase64 = mu * interv;
//    //floor to the nearest multiples of interv
//    return phase64 >> 32;
    auto interv = static_cast<int32_t>(TORUS_Q / Msize);
    int32_t mod = intModP(mu, static_cast<int32_t>(Msize));
    return mod * interv;
}

int32_t modSwitchFromTorus32(Torus in, uint32_t newMod) {
//    uint64_t interv = ((UINT64_C(1) << 63) / newMod) * 2; // width of each interval
//    uint64_t half_interval = interv / 2; // begin of the first intervall
//    uint64_t phase64 = (uint64_t(in) << 32) + half_interval;
//    //floor to the nearest multiples of interv
//    return (in >= 0) ? (phase64 / interv) : (phase64 / interv - newMod);
    auto interv = static_cast<int32_t>(TORUS_Q / newMod);
    double div = (double)in / interv;
    auto real = static_cast<int32_t>(round(div));
    return intModP(real, static_cast<int32_t>(newMod));
}

uint32_t modSwitchFromTorus32Pos(Torus in, uint32_t newMod) {
//    uint64_t interv = ((UINT64_C(1) << 63) / newMod) * 2; // width of each interval
//    uint64_t half_interval = interv / 2; // begin of the first intervall
//    uint64_t phase64 = (uint64_t(in) << 32) + half_interval;
//    //floor to the nearest multiples of interv
//    return phase64 / interv;
    auto tmp = modSwitchFromTorus32(in, newMod); // [-newMod/2, newMod/2)
    return tmp + newMod / 2;
}

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs) {
    for (auto& coeff : coeffs) {
        coeff = uniformTorusDistrib(rng);
    }
}

void initCoeffsWithGaussianNoiseSingleSample(std::vector<Torus>& coeffs, const Torus msg, const double sigma) {
    for (auto& coeff : coeffs) {
        coeff = addGaussianNoise(msg, sigma);
    }
}

void initCoeffsWithGaussianNoiseMultiSample(std::vector<Torus>& coeffs, const std::vector<Torus>& msg, const double sigma) {
    for (auto i = 0; i < coeffs.size(); i++) {
        coeffs[i] = addGaussianNoise(msg[i], sigma);
    }
}

// (x^y) % mod
Integer modPow(Integer x, Integer y, Integer mod) {
    Integer res = 1;
    while(y != 0) {
        if ((y & 1) != 0) {
            res = (Integer)(((long long)res * x) % mod);
        }
        x = (Integer)(((long long)x * x) % mod);
        y >>= 1;
    }
    return res;
}

// Calculate the modular multiplicative inverse of 'a' modulo 'mod'
Integer reciprocal(Integer a, Integer mod) {
    Integer m0 = mod, t, q;
    Integer x0 = 0, x1 = 1;

    if (mod == 1) {
        return 0;
    }

    while (a > 1) {
        q = a / mod; // q is quotient
        t = mod;
        mod = a % mod; // m is remainder now, process same as Euclid's algorithm
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

    // Make x1 positive
    if (x1 < 0) {
        x1 += m0;
    }

    return x1;
}

bool isPrime(int num) {
    if (num <= 1) {
        return false;
    }
    for (int i = 2; i <= sqrt(num); i++) {
        if (num % i == 0) {
            return false;
        }
    }
    return true;
}

void uniquePrimeFactors(std::vector<int>& result, int n) {
    if (n < 1) { throw std::invalid_argument("Invalid input"); }

    for (int i = 2, end = sqrtFloor(n); i <= end; i++) {
        if (n % i == 0) {
            result.push_back(i);
            do {
                n /= i;
            } while (n % i == 0);
            end = sqrtFloor(n);
        }
    }
    if (n > 1) { result.push_back(n); }
}

// floor(sqrt(x))
int sqrtFloor(int x) {
    if (x < 0)
        throw std::invalid_argument("Invalid input");

    int y = 0;
    for (int i = 1 << 15; i != 0; i >>= 1) {
        y |= i;
        if (y > 46340 || y * y > x)
            y ^= i;
    }
    return y;
}