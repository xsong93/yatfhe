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
            tw_factor(N, 0), N(N), Q(Ntt32((MOD - 1) / N)) {};
    explicit Ntt32_TW(int N, Ntt32 val):
            tw_factor(N, val), N(N), Q(Ntt32((MOD - 1) / N)) {};
};
struct Ntt32_iTW {
    std::vector<Ntt32> itw_factor {};
    int N {};
    Ntt32 Q {};
    Ntt32_iTW() : itw_factor(), N(), Q() {};
    explicit Ntt32_iTW(int N):
            itw_factor(N, 0), N(N), Q(Ntt32((MOD - 1) / N)) {};
    explicit Ntt32_iTW(int N, Ntt32 val):
            itw_factor(N, val), N(N), Q(Ntt32((MOD - 1) / N)) {};
};
void print_myNtt(const Ntt32_TW& TW, const Ntt32_iTW& iTW);
void printNttPoly(const NttPolynomial& in);
int clog2(int N);
Ntt32 find_primitive_root(Ntt32 P);
Ntt32 POW(Ntt32 base, Ntt32 exp, Ntt32 mod);
Ntt32 modINV(Ntt32 in);
Ntt32 modADD(Ntt32 a, Ntt32 b);
Ntt32 modSUB(Ntt32 a, Ntt32 b);
Ntt32 modMULT(Ntt32 a, Ntt32 b);
void genTW(Ntt32_TW& TW);
void geniTW(Ntt32_iTW& iTW, const Ntt32_TW& TW);
void doNTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_TW& TW_param);
void doINTT32(NttPolynomial& res, const NttPolynomial& in, const Ntt32_iTW& iTW_param);
#endif //HLS_YATFHE_MYNTT_H
