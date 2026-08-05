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
using Decomp = int8_t;
#elif defined(TORUS56)
using Torus = int64_t;
using UnsignedInteger = uint64_t;
using Decomp = int8_t;
#else
#error "torus.h: TORUS undefined"
#endif

// LWE-domain torus, selected independently of the TRLWE Torus so a narrow
// TLWE modulus (qLwe) need not pay for the full Torus width.
//   default     : LweTorus = Torus  -- qLwe may use the full Torus range
//                 (e.g. up to 56-bit in TORUS56 builds); no behavior change.
//   LWE_TORUS32 : LweTorus = int32_t -- qLwe <= 2^32; halves TLWE/KSK/TLEV
//                 storage and key-switch memory bandwidth.
// Future-proofing note: a full 64-bit qLwe fits LweTorus == int64_t, but the
// intermediate products in multTglevWithConst() (tlev.cpp) and the KSK
// accumulation (keyswitching.cpp) are only 64-bit wide. Keeping qLwe under
// ~56 bits leaves headroom; a true 64-bit modulus would need 128-bit
// intermediates there.
#if defined(LWE_TORUS32)
using LweTorus = int32_t;
using UnsignedLwe = uint32_t;
#else
using LweTorus = Torus;
using UnsignedLwe = UnsignedInteger;
#endif

using Ntt14 = uint16_t;
using Ntt16 = uint32_t;
using Ntt24 = uint32_t;
using Ntt32 = uint32_t;
using Ntt64 = uint64_t;
using NttType = uint64_t;
using Integer = Torus;
using Binary = Integer;
constexpr int32_t NUM_PRIMES{4};
constexpr int32_t NUM_HIGH_PRIMES{2};
constexpr int32_t NUM_LOW_PRIMES{2};
//const Integer QD_CRT[NUM_PRIMES] {251, 19, 17, 13};

constexpr int QD_CRT[NUM_PRIMES] {251, 241, 239, 233};
constexpr int64_t Q_60 = INT64_C(1) << 60;
constexpr int64_t Q_59 = INT64_C(1) << 59;
constexpr int64_t Q_58 = INT64_C(1) << 58;
constexpr int64_t Q_57 = INT64_C(1) << 57;
constexpr int64_t Q_56 = INT64_C(1) << 56;
constexpr int64_t Q_55 = INT64_C(1) << 55;
constexpr int64_t Q_54 = INT64_C(1) << 54;
constexpr int64_t Q_53 = INT64_C(1) << 53;
constexpr int64_t Q_52 = INT64_C(1) << 52;
constexpr int64_t Q_51 = INT64_C(1) << 51;
constexpr int64_t Q_50 = INT64_C(1) << 50;
constexpr int64_t Q_49 = INT64_C(1) << 49;
constexpr int64_t Q_48 = INT64_C(1) << 48;
constexpr int64_t Q_47 = INT64_C(1) << 47;
constexpr int64_t Q_46 = INT64_C(1) << 46;
constexpr int64_t Q_45 = INT64_C(1) << 45;
constexpr int64_t Q_44 = INT64_C(1) << 44;
constexpr int64_t Q_43 = INT64_C(1) << 43;
constexpr int64_t Q_42 = INT64_C(1) << 42;
constexpr int64_t Q_41 = INT64_C(1) << 41;
constexpr int64_t Q_40 = INT64_C(1) << 40;
constexpr int64_t Q_39 = INT64_C(1) << 39;
constexpr int64_t Q_38 = INT64_C(1) << 38;
constexpr int64_t Q_37 = INT64_C(1) << 37;
constexpr int64_t Q_36 = INT64_C(1) << 36;
constexpr int64_t Q_35 = INT64_C(1) << 35;
constexpr int64_t Q_34 = INT64_C(1) << 34;
constexpr int64_t Q_33 = INT64_C(1) << 33;
constexpr int64_t Q_32 = INT64_C(1) << 32;
constexpr int64_t Q_31 = INT64_C(1) << 31;
constexpr int64_t Q_30 = INT64_C(1) << 30;
constexpr int64_t Q_29 = INT64_C(1) << 29;
constexpr int64_t Q_28 = INT64_C(1) << 28;
constexpr int64_t Q_27 = INT64_C(1) << 27;
constexpr int64_t Q_26 = INT64_C(1) << 26;
constexpr int64_t Q_25 = INT64_C(1) << 25;
constexpr int64_t Q_24 = INT64_C(1) << 24;
constexpr int64_t Q_23 = INT64_C(1) << 23;
constexpr int64_t Q_22 = INT64_C(1) << 22;
constexpr int64_t Q_21 = INT64_C(1) << 21;
constexpr int64_t Q_20 = INT64_C(1) << 20;
constexpr int64_t Q_19 = INT64_C(1) << 19;
constexpr int64_t Q_18 = INT64_C(1) << 18;
constexpr int64_t Q_17 = INT64_C(1) << 17;
constexpr int64_t Q_16 = INT64_C(1) << 16;
constexpr int64_t Q_15 = INT64_C(1) << 15;
constexpr int64_t Q_14 = INT64_C(1) << 14;
constexpr int64_t Q_13 = INT64_C(1) << 13;
constexpr int64_t Q_12 = INT64_C(1) << 12;
constexpr uint64_t Q_64P = 18446744073709550593ULL;
constexpr uint64_t Q_63P = 9223369837831520257ULL;
constexpr uint64_t Q_62P = 4611686018427322369ULL;
constexpr uint64_t Q_61P = 2305843009211596801ULL;
constexpr uint64_t Q_60P = 1152921504606830593ULL;
constexpr uint64_t Q_59P = 576460752303421441ULL;
constexpr uint64_t Q_58P = 288230376151711681ULL;
constexpr uint64_t Q_57P = 144115188075593729ULL;
constexpr uint64_t Q_56P = 72057594036879361ULL;
constexpr uint64_t Q_55P = 36028797018963841ULL;
constexpr uint64_t Q_54P = 18014398509481729ULL;
constexpr uint64_t Q_53P = 9007199254740481ULL;
constexpr uint64_t Q_52P = 4503599627366401ULL;
constexpr uint64_t Q_51P = 2251799813684737ULL;
constexpr uint64_t Q_50P = 1125899906826241ULL;
constexpr uint64_t Q_49P = 562948879679489ULL;
constexpr uint64_t Q_48P = 281474976694273ULL;
constexpr uint64_t Q_47P = 140737488355201ULL;
constexpr uint64_t Q_46P = 70368744177601ULL;
constexpr uint64_t Q_45P = 35184372088321ULL;
constexpr uint64_t Q_44P = 17592186028033ULL;
constexpr uint64_t Q_43P = 8796093021953ULL;
constexpr uint64_t Q_42P = 4393751543809ULL;
constexpr uint64_t Q_41P = 2199023255489ULL;
constexpr uint64_t Q_40P = 1095216660481ULL;
constexpr uint64_t Q_39P = 549755809793ULL;
constexpr uint64_t Q_38P = 274877905921ULL;
constexpr uint64_t Q_37P = 137438822401ULL;
constexpr uint64_t Q_36P = 68718428161ULL;
constexpr uint64_t Q_35P = 34359214081ULL;
constexpr uint64_t Q_34P = 17175674881ULL;
constexpr uint64_t Q_33P = 8588886017ULL;
constexpr uint32_t Q_32P = 4293918721U;
constexpr uint32_t Q_31P = 2147483137U;
constexpr uint32_t Q_30P = 1073479681U;
constexpr uint32_t Q_29P = 536870849U;
constexpr uint32_t Q_28P = 268369921U;
constexpr uint32_t Q_27P = 134215681U;
constexpr uint32_t Q_26P = 67104769U;
constexpr uint32_t Q_25P = 33550337U;
constexpr uint32_t Q_24P = 16777153U;
constexpr uint32_t Q_23P = 8380417U;
constexpr uint32_t Q_22P = 4190209U;
constexpr uint32_t Q_21P = 2097153U;
constexpr uint32_t Q_20P = 1048577U;
constexpr uint32_t Q_19P = 523777U;
constexpr uint32_t Q_18P = 262145U;
constexpr uint32_t Q_17P = 131009U;
constexpr uint32_t Q_16P = 64513U;
constexpr uint32_t Q_15P = 32257U;
constexpr uint32_t Q_14P = 15361U;
constexpr uint32_t Q_13P = 7937U;
constexpr uint32_t Q_12P = 2049U;
// NTT prime for the 32-bit torus: 131070 * 2^32 + 1. Two properties earn it a
// dedicated constant instead of one of the Q_xxP above:
//   * < 2^50, so HEXL dispatches the 52-bit AVX512-IFMA NTT and the float
//     EltwiseMultMod instead of the 64-bit paths (~1.65x on the PBS phase).
//   * = 1 (mod 2^32), so an INTT wrap costs 1 unit of torus noise rather than
//     |centred(qNtt mod q)| units -- see nttWrapDelta() in yatfhe_parameters.h.
// Being a multiple of 2^32 plus one, it is also 1 (mod 2N) for every N <= 2^31,
// so the negacyclic NTT exists for any ring degree in use here.
// (Q_50P also clears both bars, with the same wrap cost as Q_60P: 16383.
//  Q_49P does NOT -- its residue mod 2^32 is -(2^30 - 1), which is fatal.)
constexpr uint64_t Q_49P_T32 = 562941363486721ULL;
constexpr int64_t Q_CRT = static_cast<int64_t>(QD_CRT[0]) * QD_CRT[1] * QD_CRT[2] * QD_CRT[3];
constexpr uint64_t BARRETT_CONSTANT = UINT64_MAX / static_cast<uint64_t>(Q_CRT);  // μ = floor(2^64 / TORUS_Q)

extern int64_t TORUS_Q;
extern int64_t LWE_Q;
extern Integer MESSAGE_P;
extern Integer INT_MAX_VALUE;
extern Integer INT_MIN_VALUE;
extern Torus TORUS_MAX;
extern Torus TORUS_MIN;
// TORUS_Q is a power of two for the PBS parameter sets (Q_56, Q_32, ...), but NOT
// in general: the CRT paths run with q = Q_CRT, a product of odd primes. When it
// is a power of two the centered residue is a sign-extension by TORUS_SHIFT, which
// lets the hot loops use longModPow2 and vectorise; otherwise they must fall back
// to longModP. Always branch on TORUS_IS_POW2 *outside* the loop -- TORUS_SHIFT is
// meaningless (and shifting by it is UB) when the flag is false.
extern bool TORUS_IS_POW2;
extern int TORUS_SHIFT;
extern Torus LWE_MAX;
extern Torus LWE_MIN;
extern NttType NTT_MAX;
extern NttType NTT_MIN;

#endif //HLS_YATFHE_TORUS_H
