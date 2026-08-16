//
// Created by Xintong Song on 2025/3/4.
//

#ifndef YATFHE_NTT32_H
#define YATFHE_NTT32_H

#include <gmp.h>
#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"

namespace NttNative32 {

    constexpr uint32_t MOD = Q_32P;
    constexpr uint32_t HALF_MOD = (MOD + 1) >> 1;
    constexpr int32_t PRIM_ROOT = 7;

    struct TwRom {
        int N{};
        std::vector<Ntt32> w_rom{};
        std::vector<Ntt32> inv_w_rom{};
        std::vector<Ntt32> phi_rom{};
        std::vector<Ntt32> inv_phi_rom{};

        TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};

        explicit TwRom(int n) :
                N(n), w_rom(n >> 1), phi_rom(n), inv_w_rom(n >> 1), inv_phi_rom(n) {};

        static void initTwRom(TwRom &twRom, const int n) {
            twRom = TwRom(n);
        }
    };

    struct TwParam {
        std::vector<std::vector<Ntt32>> tw_factor{};

        TwParam() : tw_factor() {};

        explicit TwParam(int n) :
                tw_factor(n, std::vector<Ntt32>()) {};

        static void initTwParam(TwParam &twParam, const int n) {
            twParam = TwParam(n);
        }
    };

    extern TwParam NWC_TW;
    extern TwParam NWC_ITW;
    extern TwRom TW_ROM;


    void genTW_ROM(TwRom &tw_rom);

    void genNWCparam(TwParam &nwc_tw, int n, const TwRom &tw_rom, const std::string &str);

    void applyNtt(Ntt32Polynomial &RES, const TorusPolynomial &IN);

    void applyIntt(TorusPolynomial &out, const Ntt32Polynomial &in);

    Ntt32 POW(Ntt32 BASE, Ntt32 EXP);

    Ntt32 modINV(Ntt32 in);

    Ntt32 modADD(Ntt32 a, Ntt32 b);

    Ntt32 modADDscale(Ntt32 a, Ntt32 b);

    Ntt32 modSUBscale(Ntt32 a, Ntt32 b);

    Ntt32 modSUB(Ntt32 a, Ntt32 b);

    Ntt32 modMULT(Ntt32 a, Ntt32 b);

    void initGlobalParamsNtt(int N);

    template<typename T, typename R>
    void applyNttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyNtt(out.a[row], in.a[row]);
        }
        applyNtt(out.b, in.b);
    }

    template<typename T, typename R>
    void applyInttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyIntt(out.a[row], in.a[row]);
        }
        applyIntt(out.b, in.b);
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

    template<typename RgswDftType, typename RgswType>
    void applyInttForRgsw(RgswType &out, RgswDftType &in) {
        auto level = out.l;
        auto k = out.k;
        for (size_t l = 0; l < level; l++) {
            for (size_t k1 = 0; k1 < k + 1; k1++) {
                auto &dftIn = in.trlweDftSamples[l][k1];
                auto &rgswOut = out.trlweSamples[l][k1];
                applyInttForAB(rgswOut, dftIn);
            }
        }
    }

    void calModularInnerProductNtt(Ntt32Polynomial &out, const Ntt32Polynomial &in1, const Ntt32Polynomial &in2);
}

#endif //YATFHE_NTT32_H
