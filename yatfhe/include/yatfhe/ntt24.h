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

namespace NttNative24 {

    constexpr int32_t MOD = (1 << 24) - (1 << 14) + 1;
    constexpr uint32_t HALF_MOD = (MOD + 1) >> 1;
    constexpr int32_t PRIM_ROOT = 7;
    constexpr int32_t NTT_MASK = 0xff;

    struct TwRom {
        int N{};
        std::vector<Ntt24> w_rom{};
        std::vector<Ntt24> inv_w_rom{};
        std::vector<Ntt24> phi_rom{};
        std::vector<Ntt24> inv_phi_rom{};

        TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};

        explicit TwRom(int n) :
                N(n), w_rom(n >> 1), phi_rom(n), inv_w_rom(n >> 1), inv_phi_rom(n) {};

        static void initTwRom(TwRom &twRom, const int n) {
            twRom = TwRom(n);
        }
    };

    struct TwParam {
        std::vector<std::vector<Ntt24>> tw_factor{};

        TwParam() : tw_factor() {};

        explicit TwParam(int n) :
                tw_factor(n, std::vector<Ntt24>()) {};

        static void initTwParam(TwParam &twParam, const int n) {
            twParam = TwParam(n);
        }
    };

    extern TwParam NWC_TW;
    extern TwParam NWC_ITW;
    extern TwRom TW_ROM;

//----------------------------------------------------------------------------------
    void genTW_ROM(TwRom &tw_rom);

    void genNWCparam(TwParam &nwc_tw, int n, const TwRom &tw_rom, const std::string &str);

    void applyNtt(Ntt24Polynomial &RES, const Int8Polynomial &IN);

    void applyIntt(Int8Polynomial &out, const Ntt24Polynomial &in, int q);

    Ntt24 POW(Ntt24 BASE, Ntt24 EXP);

    Ntt24 modINV(Ntt24 in);

    Ntt24 modADD(Ntt24 a, Ntt24 b);

    Ntt24 modADDscale(Ntt24 a, Ntt24 b);

    Ntt24 modSUBscale(Ntt24 a, Ntt24 b);

    Ntt24 modSUB(Ntt24 a, Ntt24 b);

    Ntt24 modMULT(Ntt24 a, Ntt24 b);

    void initGlobalParamsNtt(int N);

    template<typename T, typename R>
    void applyNttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyNtt(out.a[row], in.a[row]);
        }
        applyNtt(out.b, in.b);
    }

    template<typename T, typename R, typename U>
    void applyInttForAB(T &out, R &in, U q) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyIntt(out.a[row], in.a[row], q);
        }
        applyIntt(out.b, in.b, q);
    }

    template<typename RgswDftType, typename RgswType>
    void applyNttForRgsw(RgswDftType &out, RgswType &in) {
        auto level = in.l;
        auto k = in.k;
        for (size_t l = 0; l < level; l++) {
            for (size_t k1 = 0; k1 < k + 1; k1++) {
                auto &nttOut = out.trlweDftSamples[l][k1];
                auto &nttIn = in.trlweSamples[l][k1];
                applyNttForAB(nttOut, nttIn);
            }
        }
    }

    template<typename RgswDftType, typename RgswType, typename U>
    void applyInttForRgsw(RgswType &out, RgswDftType &in, U q) {
        auto level = out.l;
        auto k = out.k;
        for (size_t l = 0; l < level; l++) {
            for (size_t k1 = 0; k1 < k + 1; k1++) {
                auto &dftIn = in.trlweDftSamples[l][k1];
                auto &rgswOut = out.trlweSamples[l][k1];
                applyInttForAB(rgswOut, dftIn, q);
            }
        }
    }

    void calModularInnerProductNtt(Ntt24Polynomial &out, const Ntt24Polynomial &in1, const Ntt24Polynomial &in2);
}

#endif //HLS_YATFHE_NTT24_H
