//
// Created by ic on 24-4-10.
//

#ifndef HLS_YATFHE_MYNTT_H
#define HLS_YATFHE_MYNTT_H

#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "numeric_functions.h"
#include "yatfhe/torus.h"
#include <vector>
#include <gmp.h>

using Ntt32 = uint32_t;
constexpr Ntt32 MOD = 0xfff00001;
constexpr Ntt32 PRIM_ROOT = 19;
struct Ntt32_TW {
    std::vector<Ntt32> tw_factor {};
    int N {};
    Ntt32 Q {};
    Ntt32_TW() : tw_factor(), N(), Q() {};
    explicit Ntt32_TW(int N) :
            tw_factor(N, 0), N(N), Q(Ntt32((MOD - 1) / (N<<1))) {};
    explicit Ntt32_TW(int N, Ntt32 val):
            tw_factor(N, val), N(N), Q(Ntt32((MOD - 1) / (N<<1))) {};
};
struct Ntt32_PARAM {
    std::vector<Ntt32> tw_factor {};
    std::vector<Ntt32> phi_factor {};
    int N {};
    int tw_N {};
    int phi_N {};

    Ntt32_PARAM() : tw_factor(), phi_factor(), N(), tw_N(), phi_N() {};
    explicit Ntt32_PARAM(int n, int tw_n, int phi_n) :
//            tw_factor(tw_N, 0), phi_factor(phi_N, 0), tw_N(N>>1), phi_N(N), N(n){};
            tw_factor(tw_N, 0), phi_factor(phi_N, 0), tw_N(tw_n), phi_N(phi_n), N(n){};

};

struct INtt32_PARAM {
    std::vector<Ntt32> inv_tw_factor {};
    std::vector<Ntt32> inv_phi_factor {};
    int N {};
    int inv_tw_N {};
    int inv_phi_N {};

    INtt32_PARAM() : inv_tw_factor(), inv_phi_factor(), N(), inv_tw_N(), inv_phi_N() {};
    explicit INtt32_PARAM(int n, int inv_tw_n, int inv_phi_n) :
            inv_tw_factor(inv_tw_N, 0), inv_phi_factor(inv_phi_N, 0), inv_tw_N(inv_tw_n), inv_phi_N(inv_phi_n), N(n) {};
};

struct ROM {
    Ntt32_PARAM ntt_rom;
    INtt32_PARAM intt_rom;
    int N {};
    explicit ROM(int n) :
       N(n), ntt_rom(n, n, n*2), intt_rom(n, n, n*2) {};
};

struct Ntt32_iTW {
    std::vector<Ntt32> itw_factor {};
    int N {};
    Ntt32 Q {};
    Ntt32_iTW() : itw_factor(), N(), Q() {};
    explicit Ntt32_iTW(int N):
            itw_factor(N, 0), N(N), Q(Ntt32((MOD - 1) / (N<<1))) {};
    explicit Ntt32_iTW(int N, Ntt32 val):
            itw_factor(N, val), N(N), Q(Ntt32((MOD - 1) / (N<<1))) {};
};

//--------------------------------------------------------------------------------

void NWC_NTT32(NttPolynomial& RES, const NttPolynomial& IN, const Ntt32_PARAM& ntt_param);
void NWC_INTT32(NttPolynomial& RES, const NttPolynomial& IN, const INtt32_PARAM& intt_param);
void genNTT32_PARAM(Ntt32_PARAM& ntt_param);
void genINTT32_PARAM(INtt32_PARAM& intt_param, const Ntt32_PARAM& ntt_param);
void genROM(ROM& rom);
void pre_process(NttPolynomial& in, const Ntt32_PARAM& para);

//--------------------------------------------------------------------------------

void print_myNtt(const Ntt32_TW& TW, const Ntt32_iTW& iTW);
void printNttPoly(const NttPolynomial& in);
int clog2(int N);
Ntt32 find_primitive_root(Ntt32 P);
Ntt32 POW(Ntt32 base, Ntt32 exp, Ntt32 mod);
Ntt32 modINV(Ntt32 in);
Ntt32 modADD(Ntt32 a, Ntt32 b);
Ntt32 modADDscale(Ntt32 a, Ntt32 b);
Ntt32 modSUBscale(Ntt32 a, Ntt32 b);
Ntt32 modSUB(Ntt32 a, Ntt32 b);
Ntt32 modMULT(Ntt32 a, Ntt32 b);
void genTW(Ntt32_TW& TW);
void geniTW(Ntt32_iTW& iTW, const Ntt32_TW& TW);
void doNTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_TW& TW_param);
void doINTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_iTW& iTW_param);



void printPARAM(Ntt32_PARAM& ntt_param, INtt32_PARAM& intt_param);
void printROM(ROM& rom);
#endif //HLS_YATFHE_MYNTT_H
