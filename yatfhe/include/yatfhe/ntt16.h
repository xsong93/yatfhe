//
// Created by Xintong Song on 2024/5/28.
//

#ifndef HLS_YATFHE_NTT16_H
#define HLS_YATFHE_NTT16_H

#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"

constexpr int32_t MOD16 = 65537;
constexpr uint16_t HALF_MOD16 = (MOD16 + 1) >> 1;
constexpr unsigned char NTT16_MASK = 0xff;
constexpr char PRIM_ROOT16 = 3;
constexpr unsigned char POLY_MAX8 = 1 << 7;


struct TwRom16 {
    int N {};
    std::vector<Ntt16> w_rom {};
    std::vector<Ntt16> inv_w_rom {};
    std::vector<Ntt16> phi_rom {};
    std::vector<Ntt16> inv_phi_rom {};

    TwRom16() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TwRom16(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
    static void initTwRom(TwRom16& twRom, const int n) {
        twRom = TwRom16(n);
    }
};

struct TwParam16 {
    std::vector<std::vector<Ntt16>> tw_factor {};
    TwParam16(): tw_factor() {};

    explicit TwParam16(int n):
            tw_factor(n, std::vector<Ntt16>()) {};

    static void initTwParam(TwParam16& twParam, const int n) {
        twParam = TwParam16(n);
    }
};

struct Ntt16Polynomial {
    std::vector<Ntt16> coeffs {}; // N
    int N {};

    explicit Ntt16Polynomial(int N) :
            N(N),
            coeffs(N, 0) {};

    Ntt16Polynomial(int N, Ntt16 value) :
            N(N),
            coeffs(N, value) {};
};


extern TwParam16 NWC_TW16;
extern TwParam16 NWC_ITW16;
extern TwRom16 TW_ROM16;

//----------------------------------------------------------------------------------
void genTW_ROM16(TwRom16& tw_rom);
void genNWCparam16(TwParam16& nwc_tw,const int n, const TwRom16& tw_rom, const std::string& str);
void applyNtt16(Ntt16Polynomial& RES, const IntPolynomial& IN);
void applyIntt16(IntPolynomial & RES, const Ntt16Polynomial& IN);
Ntt16 POW16(Ntt16 base, Ntt16 exp, Ntt16 mod);
Ntt16 modINV16(Ntt16 in);
Ntt16 modADD16(Ntt16 a, Ntt16 b);
Ntt16 modADDscale16(Ntt16 a, Ntt16 b);
Ntt16 modSUBscale16(Ntt16 a, Ntt16 b);
Ntt16 modSUB16(Ntt16 a, Ntt16 b);
Ntt16 modMULT16(Ntt16 a, Ntt16 b);
void initGlobalParamsNtt16(int N);

#endif //HLS_YATFHE_NTT16_H
