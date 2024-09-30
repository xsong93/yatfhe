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
using Ntt14 = uint16_t;
using Ntt16 = uint32_t;
using Ntt24 = uint32_t;
using Ntt64 = uint64_t;
using UnsignedInteger = uint32_t;
using NttType = uint64_t;
using Integer = int32_t;
using Binary = Integer;
const Integer QD_CRT[4] {251, 241, 239, 233};
//const uint64_t TORUS_Q = INT64_C(1) << 32;
const int64_t TORUS_Q {static_cast<int64_t>(QD_CRT[0]) * QD_CRT[1] * QD_CRT[2] * QD_CRT[3]};
const int64_t POS_HALF_TORUS_Q = (TORUS_Q -1)>>1;
const int64_t NEG_HALF_TORUS_Q = -((TORUS_Q - 1)>>1);
//const int64_t TORUS_Q_RANGE = (TORUS_Q + 1);
const Integer INT_MAX_CRT = static_cast<Integer>(TORUS_Q / 2);
const Integer INT_MIN_CRT = -INT_MAX_CRT;
const Integer TORUS_MAX = INT_MAX_CRT;
const Integer TORUS_MIN = INT_MIN_CRT;
const uint64_t BARRETT_CONSTANT = UINT64_MAX / (uint64_t)TORUS_Q;  // μ = floor(2^64 / TORUS_Q)

#else
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Integer = int64_t;
using Binary = Integer;
const Integer TorusMax = INT64_MAX;
const Integer TorusMin = INT64_MIN;
#endif

#endif //HLS_YATFHE_TORUS_H
