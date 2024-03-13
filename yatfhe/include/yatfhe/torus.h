//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <vector>
//#include <cstdint>

using Binary = int16_t;

//#ifdef TORUS32
using Torus = int32_t;
using UnsignedInteger = uint32_t;
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
using Integer = Torus;

#endif //HLS_YATFHE_TORUS_H
