//
// Created by ic on 24-4-29.
//

#ifndef HLS_YATFHE_NTT64_H
#define HLS_YATFHE_NTT64_H

#include <vector>
#include <string>
#include <gmp.h>
#include "yatfhe/torus.h"
#include "yatfhe/ntt.h"
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"

constexpr Ntt64 MOD64 = 0xffffffff00000001;
constexpr Ntt64 HALF_MOD64 = (MOD64 + 1) >> 1;
constexpr uint32_t NTT64_MASK = 0xffffffff;
constexpr Ntt64 PRIM_ROOT64 = 7;
constexpr uint32_t POLY_MAX = 1 << 31;
const string STR_NTT = "NWC-DIT-NR-NNT";
const string STR_INTT = "NWC-DIF-RN-INNT";

struct TW_ROM {
    int N {};
    std::vector<Ntt64> w_rom {};
    std::vector<Ntt64> inv_w_rom {};
    std::vector<Ntt64> phi_rom {};
    std::vector<Ntt64> inv_phi_rom {};
    TW_ROM() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TW_ROM(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
};

struct TW_PARAM {
    std::vector<std::vector<Ntt64>> tw_factor {};
    TW_PARAM(): tw_factor() {};

    explicit TW_PARAM(int n):
            tw_factor(n, std::vector<Ntt64>()) {};

//    TW_PARAM(int n, bool init):
//    tw_factor(n, std::vector<Ntt64>()) {
//        if (init) {
//            cout<<endl;
//        }
//    };
};


//extern TW_PARAM twParam;
//extern TW_ROM twRom;

//----------------------------------------------------------------------------------
void genTW_ROM(TW_ROM& tw_rom);
void genNWCparam(TW_PARAM& nwc_tw,const int n, const TW_ROM& tw_rom, const std::string str);
void doNTT(Ntt64Polynomial& RES, const IntPolynomial & IN, const TW_PARAM& ntt_param);
void DIT_NR(Ntt64Polynomial& RES, const Ntt64Polynomial & IN, const TW_PARAM& ntt_param);
void doINTT(IntPolynomial & RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param);
void DIF_RN(Ntt64Polynomial & RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param);
int clog2(int N);
Ntt64 POW(Ntt64 base, Ntt64 exp, Ntt64 mod);
Ntt64 modINV(Ntt64 in);
Ntt64 modADD(Ntt64 a, Ntt64 b);
Ntt64 modADDscale(Ntt64 a, Ntt64 b);
Ntt64 modSUBscale(Ntt64 a, Ntt64 b);
Ntt64 modSUB(Ntt64 a, Ntt64 b);
Ntt64 modMULT(Ntt64 a, Ntt64 b);
Ntt64 modmul(Ntt64 x, Ntt64 y);
void bit_rev(std::vector<Ntt64>& x);
#endif //HLS_YATFHE_NTT64_H
