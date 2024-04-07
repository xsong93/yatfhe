//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <vector>
#include <cstdint>

//#ifdef TORUS32
using UnsignedInteger = uint32_t;
using Integer = int32_t;
using Binary = Integer;
using Torus = int32_t;
const Integer TorusMax = INT32_MAX;
const Integer TorusMin = INT32_MIN;
//#undef USE_COMPRESSED_TRLWE
//#else
//using Torus = uint64_t;
//using UnsignedInteger = uint64_t;
//#endif

/* Polynomials */
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
