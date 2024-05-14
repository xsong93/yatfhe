//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <vector>
#include <cstdint>

#define TORUS32

#ifdef TORUS32
using Torus = int32_t; // use 32-bit int as torus to handle modular arithmetic naturally
using Ntt32 = uint32_t;
using Ntt64 = uint64_t;
using UnsignedInteger = uint32_t;
using NttType = uint64_t;
using Integer = int32_t;
using Binary = Integer;
const Integer IntMax = INT32_MAX;
const Integer IntMin = INT32_MIN;
const Integer TorusMax = IntMax;
const Integer TorusMin = IntMin;
#else
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Integer = int64_t;
using Binary = Integer;
const Integer TorusMax = INT64_MAX;
const Integer TorusMin = INT64_MIN;
#endif

struct TorusPolynomial {
    std::vector<Torus> coeffs {}; // N
    int N {};

    explicit TorusPolynomial(int N) :
        coeffs(N, 0),
        N(N) {};

    TorusPolynomial(int N, Torus value) :
            coeffs(N, value),
            N(N) {};
};
struct NttPolynomial {
    std::vector<Ntt32> coeffs {}; // N
    int N {};

    NttPolynomial() : coeffs(), N() {};

    explicit NttPolynomial(int N) :
            coeffs(N, 0),
            N(N) {};

    explicit NttPolynomial(int N, Ntt32 value) :
            coeffs(N, value),
            N(N) {};
};
struct Ntt64Polynomial {
    std::vector<NttType> coeffs {};
    int N {};
    Ntt64Polynomial() : N(), coeffs() {};
    explicit Ntt64Polynomial (int n):
        coeffs(n,0), N(n) {};
};
using IntPolynomial = TorusPolynomial;
using BinPolynomial = IntPolynomial;

#endif //HLS_YATFHE_TORUS_H
