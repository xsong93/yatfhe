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
using Torus = int32_t; // use 32-bit int as torus to handle overflow naturally
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

    TorusPolynomial() : coeffs(), N() {};

    explicit TorusPolynomial(int N) :
        coeffs(N, 0),
        N(N) {};

    explicit TorusPolynomial(int N, Torus value) :
            coeffs(N, value),
            N(N) {};
};

using IntPolynomial = TorusPolynomial;
using BinPolynomial = IntPolynomial;

#endif //HLS_YATFHE_TORUS_H
