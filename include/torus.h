//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <vector>
//#include <cstdint>

//#ifdef TORUS32
using Torus = int32_t;
//#undef USE_COMPRESSED_TRLWE
//#else
//using Torus = uint64_t;
//#endif

/* Polynomials */
struct TorusPolynomial {
//    Torus* coeffs{ new Torus[N] }; // N
    std::vector<Torus> coeffs; // N
    int N{};
};

using Binary = int16_t;
using IntPolynomial = TorusPolynomial;
using Integer = Torus;

#endif //HLS_YATFHE_TORUS_H
