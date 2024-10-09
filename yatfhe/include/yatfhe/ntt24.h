//
// Created by Xintong Song on 2024/9/13.
//

#ifndef HLS_YATFHE_NTT24_H
#define HLS_YATFHE_NTT24_H

#include <gmp.h>
#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"

constexpr int32_t MOD24 = (1 << 24) - (1 << 14) + 1;
constexpr uint32_t HALF_MOD24 = (MOD24 + 1) >> 1;
constexpr int32_t PRIM_ROOT24 = 7;
constexpr int32_t NTT24_MASK = 0xff;

struct TwRom24 {
    int N {};
    std::vector<Ntt24> w_rom {};
    std::vector<Ntt24> inv_w_rom {};
    std::vector<Ntt24> phi_rom {};
    std::vector<Ntt24> inv_phi_rom {};

    TwRom24() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TwRom24(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
    static void initTwRom(TwRom24& twRom, const int n) {
        twRom = TwRom24(n);
    }
};

struct TwParam24 {
    std::vector<std::vector<Ntt24>> tw_factor {};
    TwParam24(): tw_factor() {};

    explicit TwParam24(int n):
            tw_factor(n, std::vector<Ntt24>()) {};

    static void initTwParam(TwParam24& twParam, const int n) {
        twParam = TwParam24(n);
    }
};

extern TwParam24 NWC_TW24;
extern TwParam24 NWC_ITW24;
extern TwRom24 TW_ROM24;

//----------------------------------------------------------------------------------
void genTW_ROM24(TwRom24& tw_rom);
void genNWCparam24(TwParam24& nwc_tw, int n, const TwRom24& tw_rom, const std::string& str);
void applyNtt24(Ntt24Polynomial& RES, const Int8Polynomial& IN);
void applyIntt24(Int8Polynomial& out, const Ntt24Polynomial& in, int q);
Ntt24 POW24(Ntt24 BASE, Ntt24 EXP);
Ntt24 modINV24(Ntt24 in);
Ntt24 modADD24(Ntt24 a, Ntt24 b);
Ntt24 modADDscale24(Ntt24 a, Ntt24 b);
Ntt24 modSUBscale24(Ntt24 a, Ntt24 b);
Ntt24 modSUB24(Ntt24 a, Ntt24 b);
Ntt24 modMULT24(Ntt24 a, Ntt24 b);
void initGlobalParamsNtt24(int N);

template <typename T, typename R>
void applyNttForAB24(T& out, R& in) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyNtt24(out.a[row], in.a[row]);
    }
    applyNtt24(out.b, in.b);
}

template <typename T, typename R, typename U>
void applyInttForAB24(T& out, R& in, U q) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyIntt24(out.a[row], in.a[row], q);
    }
    applyIntt24(out.b, in.b, q);
}

void calModularInnerProductNtt24(Ntt24Polynomial& out, const Ntt24Polynomial& in1, const Ntt24Polynomial& in2);

#endif //HLS_YATFHE_NTT24_H
