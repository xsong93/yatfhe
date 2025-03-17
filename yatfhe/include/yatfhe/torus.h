//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <vector>
#include <cstdint>

#define TORUS32
//#define USE_CRT

#ifdef TORUS32
using Torus = int32_t; // use 32-bit int as torus to handle modular arithmetic naturally
using Ntt14 = uint16_t;
using Ntt16 = uint32_t;
using Ntt24 = uint32_t;
using Ntt32 = uint32_t;
using Ntt64 = uint64_t;
using UnsignedInteger = uint32_t;
using NttType = uint64_t;
using Integer = int32_t;
using Binary = Integer;
const int32_t NUM_PRIMES{4};
const int32_t NUM_HIGH_PRIMES{2};
const int32_t NUM_LOW_PRIMES{2};
//const Integer QD_CRT[NUM_PRIMES] {251, 19, 17, 13};
const Integer QD_CRT[NUM_PRIMES] {251, 241, 239, 233};
const int64_t Q_32 = INT64_C(1) << 32;
const uint32_t Q_32P = 16760833;
const int64_t Q_CRT = static_cast<int64_t>(QD_CRT[0]) * QD_CRT[1] * QD_CRT[2] * QD_CRT[3];
const uint64_t BARRETT_CONSTANT = UINT64_MAX / (uint64_t)Q_CRT;  // μ = floor(2^64 / TORUS_Q)

extern int64_t TORUS_Q;
extern Integer INT_MAX_VALUE;
extern Integer INT_MIN_VALUE;
extern Integer TORUS_MAX;
extern Integer TORUS_MIN;

#else
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Integer = int64_t;
using Binary = Integer;
#endif

#endif //HLS_YATFHE_TORUS_H
