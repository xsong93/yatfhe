//
// Created by ic on 24-4-29.
//

#ifndef HLS_YATFHE_NTT64_H
#define HLS_YATFHE_NTT64_H

#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"
#include "yatfhe/torus.h"
#include <vector>
#include <string>
#include <gmp.h>

using ntt64 = uint64_t;
constexpr ntt64 MOD64 = 0xffffffff00000001;
constexpr ntt64 HALF_MOD64 = (MOD64 + 1) >> 1;
constexpr uint32_t mask =  0xffffffff;
constexpr ntt64 PRIM_ROOT64 = 7;
constexpr uint32_t poly_max = 1<<31;

struct TW_PARAM {
    std::vector<std::vector<ntt64>> tw_factor {};
    TW_PARAM(): tw_factor() {};
    explicit TW_PARAM(int n):
            tw_factor(n, std::vector<ntt64>()) {};
};
struct TW_ROM {
    int N {};
    std::vector<ntt64> w_rom {};
    std::vector<ntt64> inv_w_rom {};
    std::vector<ntt64> phi_rom {};
    std::vector<ntt64> inv_phi_rom {};
    TW_ROM() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TW_ROM(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
};

//----------------------------------------------------------------------------------
void genTW_ROM(TW_ROM& tw_rom);
void genNWCparam(TW_PARAM& nwc_tw,const int n, const TW_ROM& tw_rom, const std::string str);
void doNTT(Ntt64Polynomial& RES, const IntPolynomial & IN, const TW_PARAM& ntt_param);
void DIT_NR(Ntt64Polynomial& RES, const Ntt64Polynomial & IN, const TW_PARAM& ntt_param);
void doINTT(IntPolynomial & RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param);
void DIF_RN(Ntt64Polynomial & RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param);
int clog2(int N);
ntt64 POW(ntt64 base, ntt64 exp, ntt64 mod);
ntt64 modINV(ntt64 in);
ntt64 modADD(ntt64 a, ntt64 b);
ntt64 modADDscale(ntt64 a, ntt64 b);
ntt64 modSUBscale(ntt64 a, ntt64 b);
ntt64 modSUB(ntt64 a, ntt64 b);
ntt64 modMULT(ntt64 a, ntt64 b);
void bit_rev(std::vector<ntt64>& x);
#endif //HLS_YATFHE_NTT64_H
