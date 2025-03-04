//
// Created by Xintong Song on 2025/3/4.
//

#ifndef HLS_YATFHE_NTT32_H
#define HLS_YATFHE_NTT32_H

#include <gmp.h>
#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"

constexpr uint32_t MOD32 = Q_32P;
constexpr uint32_t HALF_MOD32 = (MOD32 + 1) >> 1;
constexpr int32_t PRIM_ROOT32 = 7;

struct TwRom32 {
    int N {};
    std::vector<Ntt32> w_rom {};
    std::vector<Ntt32> inv_w_rom {};
    std::vector<Ntt32> phi_rom {};
    std::vector<Ntt32> inv_phi_rom {};

    TwRom32() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};
    explicit TwRom32(int n) :
            N(n), w_rom(n>>1), phi_rom(n), inv_w_rom(n>>1), inv_phi_rom(n) {};
    static void initTwRom(TwRom32& twRom, const int n) {
        twRom = TwRom32(n);
    }
};

struct TwParam32 {
    std::vector<std::vector<Ntt32>> tw_factor {};
    TwParam32(): tw_factor() {};

    explicit TwParam32(int n):
            tw_factor(n, std::vector<Ntt32>()) {};

    static void initTwParam(TwParam32& twParam, const int n) {
        twParam = TwParam32(n);
    }
};

extern TwParam32 NWC_TW32;
extern TwParam32 NWC_ITW32;
extern TwRom32 TW_ROM32;

//----------------------------------------------------------------------------------
void genTW_ROM32(TwRom32& tw_rom);
void genNWCparam32(TwParam32& nwc_tw, int n, const TwRom32& tw_rom, const std::string& str);
void applyNtt32(Ntt32Polynomial& RES, const TorusPolynomial& IN);
void applyIntt32(TorusPolynomial& out, const Ntt32Polynomial& in);
Ntt32 POW32(Ntt32 BASE, Ntt32 EXP);
Ntt32 modINV32(Ntt32 in);
Ntt32 modADD32(Ntt32 a, Ntt32 b);
Ntt32 modADDscale32(Ntt32 a, Ntt32 b);
Ntt32 modSUBscale32(Ntt32 a, Ntt32 b);
Ntt32 modSUB32(Ntt32 a, Ntt32 b);
Ntt32 modMULT32(Ntt32 a, Ntt32 b);
void initGlobalParamsNtt32(int N);

template <typename RgswDftType, typename RgswType>
void applyNttForRgsw32(RgswDftType& out, RgswType& in) {
    auto level = in.l;
    auto k = in.k;
    for (size_t l = 0; l < level; l++) {
        for (size_t k1 = 0; k1 < k + 1; k1++) {
            auto& nttOut = out.trlweDftSamples[l][k1];
            auto& nttIn = in.trlweSamples[l][k1];
            applyNttForAB32(nttOut, nttIn);
        }
    }
}

template <typename RgswDftType, typename RgswType>
void applyInttForRgsw32(RgswType& out, RgswDftType& in) {
    auto level = out.l;
    auto k = out.k;
    for (size_t l = 0; l < level; l++) {
        for (size_t k1 = 0; k1 < k + 1; k1++) {
            auto& dftIn = in.trlweDftSamples[l][k1];
            auto& rgswOut = out.trlweSamples[l][k1];
            applyInttForAB32(rgswOut, dftIn);
        }
    }
}

template <typename T, typename R>
void applyNttForAB32(T& out, R& in) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyNtt32(out.a[row], in.a[row]);
    }
    applyNtt32(out.b, in.b);
}

template <typename T, typename R>
void applyInttForAB32(T& out, R& in) {
    for (auto row = 0; row < in.a.size(); row++) {
        applyIntt32(out.a[row], in.a[row]);
    }
    applyIntt32(out.b, in.b);
}

void calModularInnerProductNtt32(Ntt32Polynomial& out, const Ntt32Polynomial& in1, const Ntt32Polynomial& in2);

#endif //HLS_YATFHE_NTT32_H
