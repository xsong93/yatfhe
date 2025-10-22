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
uniform_int_distribution<Integer> ternaryDistrib(-1, 1);

int calLogBase2(int N) {
    int res = 0;
    while (N >>= 1){
        res ++;
    }
    return res;
}

uniform_int_distribution<Torus>& uniformTorusDistrib(Torus min, Torus max) {
    static uniform_int_distribution<Torus> instance(min, max);
    return instance;
}

Integer genIntUniformDist(const Integer lowerBound, const Integer upperBound) {
    uniform_int_distribution<Integer> uniformIntDistrib(lowerBound, upperBound);
    return uniformIntDistrib(rng);
}

uint64_t genUInt64UniformDist(const uint64_t lowerBound, const uint64_t upperBound) {
    uniform_int_distribution<uint64_t> uniformIntDistrib(lowerBound, upperBound);
    return uniformIntDistrib(rng);
}

// Gaussian sample centered in message, with standard deviation sigma
Torus addGaussianNoise(Torus message, const double sigma, const int64_t torusQ) {
    normal_distribution normalDistribution(0.0, sigma);
    double e = normalDistribution(rng) * static_cast<double>(torusQ);
    auto err = static_cast<Torus>(e);
    Torus tmp = addTorus(torusQ, message, err);
    if ((message > 0 && tmp < 0) || (message < 0 && tmp > 0)) { // handle overflow
        return subTorus(torusQ, message, err);
    }
    return tmp;
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
    return modP / double(torusBase);
}

Torus roundTorusGeneralError(const Torus in, const int torusBase, const int64_t q) {
    auto t = modSwitchFromTorusGeneral(in, torusBase, q);
    auto t2 = modSwitchToTorusGeneral(t, torusBase, q);
    return t2;
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
    if (b > (p - 1) / 2) {
        b -= p;
    } else if (b < - p / 2) {
        b += p;
    }
    return b;
}

int64_t longModP(const int64_t a, const int64_t p) {
    auto b = a % p;
    if (b > (p - 1) / 2) {
        b -= p;
    } else if (b < - p / 2) {
        b += p;
    }
    return b;
}

// a * b = (k*a_hi + a_lo) * (k*b_hi + b_lo)
// k = 2^32
// https://stackoverflow.com/a/31662911/6553631
uint64_t mulhi64(uint64_t a, uint64_t b) {
    const uint64_t a_lo = static_cast<uint32_t>(a);
    const uint64_t a_hi = a >> 32;
    const uint64_t b_lo = static_cast<uint32_t>(b);
    const uint64_t b_hi = b >> 32;

    const uint64_t p0 = a_lo * b_lo;
    const uint64_t p1 = a_lo * b_hi;
    const uint64_t p2 = a_hi * b_lo;
    const uint64_t p3 = a_hi * b_hi;

    const uint32_t cy = ((p0 >> 32) + static_cast<uint32_t>(p1) + static_cast<uint32_t>(p2)) >> 32;

    return p3 + (p1 >> 32) + (p2 >> 32) + cy;
}

int64_t barrettReduceT32(int64_t in) {
    auto x = static_cast<uint64_t>(in < 0 ? -in : in);
    uint64_t r = mulhi64(x, BARRETT_CONSTANT);
    r = x - r * TORUS_Q;
    if (r >= TORUS_Q) {
        r -= TORUS_Q;
    }
    return in < 0 ? -static_cast<int64_t>(r) : static_cast<int64_t>(r);
}

int64_t montgomoryReduceT32(int64_t in) {
    auto x = static_cast<uint64_t>(in < 0 ? -in : in);
//    uint32_t xModR = x & 0xFFFFFFFF;
    uint64_t u = (x * 926404979) & 0xFFFFFFFF;
    uint64_t r = (x + u * TORUS_Q) >> 32;
    if (r >= TORUS_Q) {
        r -= TORUS_Q;
    }
    return in < 0 ? -static_cast<int64_t>(r) : static_cast<int64_t>(r);
}

Torus subTorus(int64_t q, Torus in1, Torus in2) {
    if (q == Q_32) {
        return in1 - in2;
    }
    auto tmp = static_cast<int64_t>(in1) - static_cast<int64_t>(in2);
//    return (tmp > TORUS_MAX) ? static_cast<Torus>(tmp - TORUS_Q) : static_cast<Torus>((tmp < TORUS_MIN) ? (TORUS_Q + tmp) : tmp);
    return (Torus)longModP(tmp, q);
}

Torus multTorus(int64_t q, Torus in1, Torus in2) {
    if (q == Q_32) {
        return in1 * in2;
    }
    return static_cast<Torus>(longModP(static_cast<int64_t>(in1) * static_cast<int64_t>(in2), q));
}

// Multiplicative inverse modulo p
int64_t modInverse(int64_t a, int64_t mod) {
    int64_t m0 = mod, t, q;
    int64_t x0 = 0, x1 = 1;

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

    return x1;
}

Torus modSwitchToTorusGeneral(int32_t mu, uint32_t mSize, int64_t torusQ) {
    auto interv = static_cast<int32_t>(torusQ / mSize);
    int32_t mod = intModP(mu, static_cast<int32_t>(mSize));
    return mod * interv;
}

int32_t modSwitchFromTorusGeneral(Torus in, uint32_t newMod, int64_t torusQ) {
    auto interv = static_cast<int32_t>(torusQ / newMod);
    double div = (double)in / interv;
    auto real = static_cast<int32_t>(round(div));
    return intModP(real, static_cast<int32_t>(newMod));
}

Torus modSwitchToTorus32(int32_t mu, uint32_t mSize) {
    auto interv = static_cast<int32_t>(TORUS_Q / mSize);
    int32_t mod = intModP(mu, static_cast<int32_t>(mSize));
    return mod * interv;
}

int32_t modSwitchFromTorus32(Torus in, uint32_t newMod) {
    auto interv = static_cast<int32_t>(TORUS_Q / newMod);
    double div = (double)in / interv;
    auto real = static_cast<int32_t>(round(div));
    return intModP(real, static_cast<int32_t>(newMod));
}

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, Torus min, Torus max) {
    for (auto& coeff : coeffs) {
        coeff = uniformTorusDistrib(min, max)(rng);
    }
}

void initCoeffsWithGaussianNoiseSingleSample(std::vector<Torus>& coeffs, const Torus msg, const double sigma,
                                             const int64_t torusQ) {
    for (auto& coeff : coeffs) {
        coeff = addGaussianNoise(msg, sigma, torusQ);
    }
}

void initCoeffsWithGaussianNoiseMultiSample(std::vector<Torus>& coeffs, const std::vector<Torus>& msg,
                                            const double sigma, const int64_t torusQ) {
    for (auto i = 0; i < coeffs.size(); i++) {
        coeffs[i] = addGaussianNoise(msg[i], sigma, torusQ);
    }
}