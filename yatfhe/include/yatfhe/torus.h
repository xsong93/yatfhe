//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TORUS_H
#define HLS_YATFHE_TORUS_H

#include <cstdlib>
#include <cstdint>

#if defined(TORUS32)
using Torus = int32_t; // use 32-bit int as torus to handle modular arithmetic naturally
using UnsignedInteger = uint32_t;
using Decomp = int16_t;
#elif defined(TORUS33)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int16_t;
#elif defined(TORUS35)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int8_t;
#elif defined(TORUS40)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int16_t;
#elif defined(TORUS42)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int8_t;
#else
#error "torus.h: TORUS undefined"
#endif

using Ntt14 = uint16_t;
using Ntt16 = uint32_t;
using Ntt24 = uint32_t;
using Ntt32 = uint32_t;
using Ntt64 = uint64_t;
using NttType = uint64_t;
using Integer = int32_t; // small "message-space" values (plaintext digits, secret-key bits, CRT primes) — independent of torus width
using Binary = Integer;
const int32_t NUM_PRIMES{4};
const int32_t NUM_HIGH_PRIMES{2};
const int32_t NUM_LOW_PRIMES{2};
//const Integer QD_CRT[NUM_PRIMES] {251, 19, 17, 13};
const Integer QD_CRT[NUM_PRIMES] {251, 241, 239, 233};
const int64_t Q_42 = INT64_C(1) << 42;
const int64_t Q_40 = INT64_C(1) << 40;
const int64_t Q_32 = INT64_C(1) << 32;
const int64_t Q_33 = INT64_C(1) << 33;
const int64_t Q_35 = INT64_C(1) << 35;
const int64_t Q_27 = INT64_C(1) << 27;
const int64_t Q_26 = INT64_C(1) << 26;
const int64_t Q_25 = INT64_C(1) << 25;
const int64_t Q_20 = INT64_C(1) << 20;
const int64_t Q_12 = INT64_C(1) << 12;
const uint32_t Q_32P = 16760833;
const uint64_t Q_50P = 1125899906826241;
const uint64_t Q_64P = 0xffffffff00000001;
const int64_t Q_CRT = static_cast<int64_t>(QD_CRT[0]) * QD_CRT[1] * QD_CRT[2] * QD_CRT[3];
const uint64_t BARRETT_CONSTANT = UINT64_MAX / (uint64_t)Q_CRT;  // μ = floor(2^64 / TORUS_Q)

extern int64_t TORUS_Q;
extern int64_t LWE_Q;
extern Integer MESSAGE_P;
extern Integer INT_MAX_VALUE;
extern Integer INT_MIN_VALUE;
extern Torus TORUS_MAX;
extern Torus TORUS_MIN;
extern Torus LWE_MAX;
extern Torus LWE_MIN;
extern NttType NTT_MAX;
extern NttType NTT_MIN;

#endif //HLS_YATFHE_TORUS_H
