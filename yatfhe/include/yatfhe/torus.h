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
#elif defined(TORUS36)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int8_t;
#elif defined(TORUS56)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int16_t;
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
const int64_t Q_60 = INT64_C(1) << 60;
const int64_t Q_59 = INT64_C(1) << 59;
const int64_t Q_58 = INT64_C(1) << 58;
const int64_t Q_57 = INT64_C(1) << 57;
const int64_t Q_56 = INT64_C(1) << 56;
const int64_t Q_55 = INT64_C(1) << 55;
const int64_t Q_54 = INT64_C(1) << 54;
const int64_t Q_53 = INT64_C(1) << 53;
const int64_t Q_52 = INT64_C(1) << 52;
const int64_t Q_51 = INT64_C(1) << 51;
const int64_t Q_50 = INT64_C(1) << 50;
const int64_t Q_49 = INT64_C(1) << 49;
const int64_t Q_48 = INT64_C(1) << 48;
const int64_t Q_47 = INT64_C(1) << 47;
const int64_t Q_46 = INT64_C(1) << 46;
const int64_t Q_45 = INT64_C(1) << 45;
const int64_t Q_44 = INT64_C(1) << 44;
const int64_t Q_43 = INT64_C(1) << 43;
const int64_t Q_42 = INT64_C(1) << 42;
const int64_t Q_41 = INT64_C(1) << 41;
const int64_t Q_40 = INT64_C(1) << 40;
const int64_t Q_39 = INT64_C(1) << 39;
const int64_t Q_38 = INT64_C(1) << 38;
const int64_t Q_37 = INT64_C(1) << 37;
const int64_t Q_36 = INT64_C(1) << 36;
const int64_t Q_35 = INT64_C(1) << 35;
const int64_t Q_34 = INT64_C(1) << 34;
const int64_t Q_33 = INT64_C(1) << 33;
const int64_t Q_32 = INT64_C(1) << 32;
const int64_t Q_31 = INT64_C(1) << 31;
const int64_t Q_30 = INT64_C(1) << 30;
const int64_t Q_29 = INT64_C(1) << 29;
const int64_t Q_28 = INT64_C(1) << 28;
const int64_t Q_27 = INT64_C(1) << 27;
const int64_t Q_26 = INT64_C(1) << 26;
const int64_t Q_25 = INT64_C(1) << 25;
const int64_t Q_24 = INT64_C(1) << 24;
const int64_t Q_23 = INT64_C(1) << 23;
const int64_t Q_22 = INT64_C(1) << 22;
const int64_t Q_21 = INT64_C(1) << 21;
const int64_t Q_20 = INT64_C(1) << 20;
const int64_t Q_19 = INT64_C(1) << 19;
const int64_t Q_18 = INT64_C(1) << 18;
const int64_t Q_17 = INT64_C(1) << 17;
const int64_t Q_16 = INT64_C(1) << 16;
const int64_t Q_15 = INT64_C(1) << 15;
const int64_t Q_14 = INT64_C(1) << 14;
const int64_t Q_13 = INT64_C(1) << 13;
const int64_t Q_12 = INT64_C(1) << 12;
const uint64_t Q_64P = 18446744073709550593ULL;
const uint64_t Q_63P = 9223369837831520257ULL;
const uint64_t Q_62P = 4611686018427322369ULL;
const uint64_t Q_61P = 2305843009211596801ULL;
const uint64_t Q_60P = 1152921504606830593ULL;
const uint64_t Q_59P = 576460752303421441ULL;
const uint64_t Q_58P = 288230376151711681ULL;
const uint64_t Q_57P = 144115188075593729ULL;
const uint64_t Q_56P = 72057594036879361ULL;
const uint64_t Q_55P = 36028797018963841ULL;
const uint64_t Q_54P = 18014398509481729ULL;
const uint64_t Q_53P = 9007199254740481ULL;
const uint64_t Q_52P = 4503599627366401ULL;
const uint64_t Q_51P = 2251799813684737ULL;
const uint64_t Q_50P = 1125899906826241ULL;
const uint64_t Q_49P = 562948879679489ULL;
const uint64_t Q_48P = 281474976694273ULL;
const uint64_t Q_47P = 140737488355201ULL;
const uint64_t Q_46P = 70368744177601ULL;
const uint64_t Q_45P = 35184372088321ULL;
const uint64_t Q_44P = 17592186028033ULL;
const uint64_t Q_43P = 8796093021953ULL;
const uint64_t Q_42P = 4393751543809ULL;
const uint64_t Q_41P = 2199023255489ULL;
const uint64_t Q_40P = 1095216660481ULL;
const uint64_t Q_39P = 549755809793ULL;
const uint64_t Q_38P = 274877905921ULL;
const uint64_t Q_37P = 137438822401ULL;
const uint64_t Q_36P = 68718428161ULL;
const uint64_t Q_35P = 34359214081ULL;
const uint64_t Q_34P = 17175674881ULL;
const uint64_t Q_33P = 8588886017ULL;
const uint32_t Q_32P = 4293918721U;
const uint32_t Q_31P = 2147483137U;
const uint32_t Q_30P = 1073479681U;
const uint32_t Q_29P = 536870849U;
const uint32_t Q_28P = 268369921U;
const uint32_t Q_27P = 134215681U;
const uint32_t Q_26P = 67104769U;
const uint32_t Q_25P = 33550337U;
const uint32_t Q_24P = 16777153U;
const uint32_t Q_23P = 8380417U;
const uint32_t Q_22P = 4190209U;
const uint32_t Q_21P = 2097153U;
const uint32_t Q_20P = 1048577U;
const uint32_t Q_19P = 523777U;
const uint32_t Q_18P = 262145U;
const uint32_t Q_17P = 131009U;
const uint32_t Q_16P = 64513U;
const uint32_t Q_15P = 32257U;
const uint32_t Q_14P = 15361U;
const uint32_t Q_13P = 7937U;
const uint32_t Q_12P = 2049U;
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
