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
const Integer IntMax = INT32_MAX;
const Integer IntMin = INT32_MIN;
const Integer TorusMax = IntMax;
const Integer TorusMin = IntMin;
const Integer QdCRT[4] {251, 241, 239, 233};
const Integer IntMaxCRT = 1684281158;
const Integer IntMinCRT = -IntMaxCRT - 1;
const Integer TorusMaxCRT = IntMaxCRT;
const Integer TorusMinCRT = IntMinCRT;
#else
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Integer = int64_t;
using Binary = Integer;
const Integer TorusMax = INT64_MAX;
const Integer TorusMin = INT64_MIN;
#endif

#endif //HLS_YATFHE_TORUS_H
