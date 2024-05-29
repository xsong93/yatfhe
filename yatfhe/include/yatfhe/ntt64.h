//
// Created by ic on 24-4-29.
//

#ifndef HLS_YATFHE_NTT64_H
#define HLS_YATFHE_NTT64_H

#include <vector>
#include <string>
#include <gmp.h>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"

constexpr Ntt64 MOD64 = 0xffffffff00000001;
constexpr Ntt64 HALF_MOD64 = (MOD64 + 1) >> 1;
constexpr uint32_t NTT64_MASK = 0xffffffff;
constexpr Ntt64 PRIM_ROOT64 = 7;
constexpr uint32_t POLY_MAX32 = 1 << 31;

struct TwRom {
    int N {};
    std::vector<Ntt64> w_rom {};
    std::vector<Ntt64> inv_w_rom {};
    std::vector<Ntt64> phi_rom {};
    std::vector<Ntt64> inv_phi_rom {};

    TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TwRom(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
    static void initTwRom(TwRom& twRom, const int n) {
        twRom = TwRom(n);
    }
};

struct TwParam {
    std::vector<std::vector<Ntt64>> tw_factor {};
    TwParam(): tw_factor() {};

    explicit TwParam(int n):
            tw_factor(n, std::vector<Ntt64>()) {};

    static void initTwParam(TwParam& twParam, const int n) {
        twParam = TwParam(n);
    }
};


extern TwParam NWC_TW64;
extern TwParam NWC_ITW64;
extern TwRom TW_ROM64;

//----------------------------------------------------------------------------------
void genTW_ROM64(TwRom& tw_rom);
void genNWCparam64(TwParam& nwc_tw,const int n, const TwRom& tw_rom, const std::string& str);
void doNTT64(Ntt64Polynomial& RES, const IntPolynomial & IN, const TwParam& ntt_param);
void doNTT64(Ntt64Polynomial& RES, const IntPolynomial& IN);
void DIT_NR64(Ntt64Polynomial& RES, const Ntt64Polynomial & IN, const TwParam& ntt_param);
void doINTT64(IntPolynomial & RES, const Ntt64Polynomial& IN, const TwParam& intt_param);
void doINTT64(IntPolynomial & RES, const Ntt64Polynomial& IN);
void DIF_RN64(Ntt64Polynomial & RES, const Ntt64Polynomial& IN, const TwParam& intt_param);
int clog2(int N);
Ntt64 POW64(Ntt64 base, Ntt64 exp, Ntt64 mod);
Ntt64 modINV64(Ntt64 in);
Ntt64 modADD64(Ntt64 a, Ntt64 b);
Ntt64 modADDscale64(Ntt64 a, Ntt64 b);
Ntt64 modSUBscale64(Ntt64 a, Ntt64 b);
Ntt64 modSUB64(Ntt64 a, Ntt64 b);
Ntt64 modMULT64(Ntt64 a, Ntt64 b);
Ntt64 modmul64(Ntt64 x, Ntt64 y);
void initGlobalParamsNtt64(int N);

#endif //HLS_YATFHE_NTT64_H
