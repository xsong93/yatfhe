//
// Created by Xintong Song on 2024/5/28.
//

#ifndef YATFHE_NTT14_H
#define YATFHE_NTT14_H

#include <gmp.h>
#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"

namespace NttNative14 {

    constexpr int32_t MOD = 12289;
    constexpr uint16_t HALF_MOD = (MOD + 1) >> 1;
    constexpr int8_t PRIM_ROOT = 11;
    const uint8_t NTT_MASK = 0xff;

    struct TwRom {
        int N{};
        std::vector<Ntt14> w_rom{};
        std::vector<Ntt14> inv_w_rom{};
        std::vector<Ntt14> phi_rom{};
        std::vector<Ntt14> inv_phi_rom{};

        TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};

        explicit TwRom(int n) :
                N(n), w_rom(n >> 1), phi_rom(n), inv_w_rom(n >> 1), inv_phi_rom(n) {};

        static void initTwRom(TwRom &twRom, const int n) {
            twRom = TwRom(n);
        }
    };

    struct TwParam {
        std::vector<std::vector<Ntt14>> tw_factor{};

        TwParam() : tw_factor() {};

        explicit TwParam(int n) :
                tw_factor(n, std::vector<Ntt14>()) {};

        static void initTwParam(TwParam &twParam, const int n) {
            twParam = TwParam(n);
        }
    };

    extern TwParam NWC_TW14;
    extern TwParam NWC_ITW14;
    extern TwRom TW_ROM14;

    void genTW_ROM(TwRom &tw_rom);

    void genNWCparam(TwParam &nwc_tw, const int n, const TwRom &tw_rom, const std::string &str);

    void applyNtt(Ntt14Polynomial &RES, const IntPolynomial &IN);

    void applyNttPoly8(Ntt14Polynomial &RES, const Int8Polynomial &IN);

    void applyIntt(IntPolynomial &RES, const Ntt14Polynomial &IN);

    void applyInttPoly8(Int8Polynomial &RES, const Ntt14Polynomial &IN);

    Ntt14 POW(Ntt14 BASE, Ntt14 EXP);

    Ntt14 modINV(Ntt14 in);

    Ntt14 modADD(Ntt14 a, Ntt14 b);

    Ntt14 modADDscale(Ntt14 a, Ntt14 b);

    Ntt14 modSUBscale(Ntt14 a, Ntt14 b);

    Ntt14 modSUB(Ntt14 a, Ntt14 b);

    Ntt14 modMULT(Ntt14 a, Ntt14 b);

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

    void calModularInnerProductNtt(Ntt14Polynomial &b, const Ntt14Polynomial &a, const Ntt14Polynomial &s);
}

#endif //YATFHE_NTT14_H
