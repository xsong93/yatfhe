//
// Created by Xintong Song on 2024/5/28.
//

#ifndef HLS_YATFHE_NTT14_H
#define HLS_YATFHE_NTT14_H

#include <gmp.h>
#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"

constexpr int32_t MOD14 = 12289;
constexpr uint16_t HALF_MOD14 = (MOD14 + 1) >> 1;
constexpr int8_t PRIM_ROOT14 = 11;
const uint8_t NTT14_MASK = 0xff;

struct TwRom14 {
    int N {};
    std::vector<Ntt14> w_rom {};
    std::vector<Ntt14> inv_w_rom {};
    std::vector<Ntt14> phi_rom {};
    std::vector<Ntt14> inv_phi_rom {};

    TwRom14() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TwRom14(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
    static void initTwRom(TwRom14& twRom, const int n) {
        twRom = TwRom14(n);
    }
};

struct TwParam14 {
    std::vector<std::vector<Ntt14>> tw_factor {};
    TwParam14(): tw_factor() {};

    explicit TwParam14(int n):
            tw_factor(n, std::vector<Ntt14>()) {};

    static void initTwParam(TwParam14& twParam, const int n) {
        twParam = TwParam14(n);
    }
};

extern TwParam14 NWC_TW14;
extern TwParam14 NWC_ITW14;
extern TwRom14 TW_ROM14;

//----------------------------------------------------------------------------------
void genTW_ROM14(TwRom14& tw_rom);
void genNWCparam14(TwParam14& nwc_tw,const int n, const TwRom14& tw_rom, const std::string& str);
void applyNtt14(Ntt14Polynomial& RES, const IntPolynomial& IN);
void applyIntt14(IntPolynomial & RES, const Ntt14Polynomial& IN);
Ntt14 POW14(Ntt14 BASE, Ntt14 EXP);
Ntt14 modINV14(Ntt14 in);
Ntt14 modADD14(Ntt14 a, Ntt14 b);
Ntt14 modADDscale14(Ntt14 a, Ntt14 b);
Ntt14 modSUBscale14(Ntt14 a, Ntt14 b);
Ntt14 modSUB14(Ntt14 a, Ntt14 b);
Ntt14 modMULT14(Ntt14 a, Ntt14 b);
void initGlobalParamsNtt14(int N);

#endif //HLS_YATFHE_NTT14_H
