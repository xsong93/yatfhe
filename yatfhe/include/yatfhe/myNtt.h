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
    explicit Ntt32_PARAM(int n) :
            tw_factor(tw_N, 0), phi_factor(phi_N, 0), tw_N(N>>1), phi_N(N), N(n){};
    explicit Ntt32_PARAM(int n, Ntt32 val):
            tw_factor(tw_N, val), phi_factor(phi_N, val), tw_N(N>>1), phi_N(N), N(n) {};
};

struct INtt32_PARAM {
    std::vector<Ntt32> itw_factor {};
    std::vector<Ntt32> iphi_factor {};
    int N {};
    int itw_N {};
    int iphi_N {};

    INtt32_PARAM() : itw_factor(), iphi_factor(), N(), itw_N(), iphi_N() {};
    explicit INtt32_PARAM(int n) :
            itw_factor(itw_N, 0), iphi_factor(iphi_N, 0), itw_N(N>>1), iphi_N(N), N(n) {};
    explicit INtt32_PARAM(int n, Ntt32 val):
            itw_factor(itw_N, val), iphi_factor(iphi_N, val), itw_N(N>>1), iphi_N(N), N(n) {};
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

void NWC_NTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_PARAM& ntt_param);
void NWC_INTT32(NttPolynomial& res, const NttPolynomial& in, const INtt32_PARAM& intt_param);
void genNTT32_PARAM(Ntt32_PARAM& ntt_param);
void genINTT32_PARAM(INtt32_PARAM& intt_param, const Ntt32_PARAM& ntt_param);

//--------------------------------------------------------------------------------

void print_myNtt(const Ntt32_TW& TW, const Ntt32_iTW& iTW);
void printNttPoly(const NttPolynomial& in);
int clog2(int N);
Ntt32 find_primitive_root(Ntt32 P);
Ntt32 POW(Ntt32 base, Ntt32 exp, Ntt32 mod);
Ntt32 modINV(Ntt32 in);
Ntt32 modADD(Ntt32 a, Ntt32 b);
Ntt32 modADDscale(Ntt32 a, Ntt32 b, bool isINTT);
Ntt32 modSUBscale(Ntt32 a, Ntt32 b, bool isINTT);
Ntt32 modSUB(Ntt32 a, Ntt32 b);
Ntt32 modMULT(Ntt32 a, Ntt32 b);
void genTW(Ntt32_TW& TW);
void geniTW(Ntt32_iTW& iTW, const Ntt32_TW& TW);
void doNTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_TW& TW_param);
void doINTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_iTW& iTW_param);
#endif //HLS_YATFHE_MYNTT_H
