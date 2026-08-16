//
// Created by Xintong Song on 2024/5/28.
//

#ifndef YATFHE_NTT16_H
#define YATFHE_NTT16_H

#include <gmp.h>
#include <vector>
#include <string>
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"

namespace NttNative16 {

    constexpr int32_t MOD = 65537;
    constexpr uint16_t HALF_MOD = (MOD + 1) >> 1;
    constexpr uint8_t NTT_MASK = 0xff;
    constexpr int8_t PRIM_ROOT = 3;

    struct TwRom {
        int N{};
        std::vector<Ntt16> w_rom{};
        std::vector<Ntt16> inv_w_rom{};
        std::vector<Ntt16> phi_rom{};
        std::vector<Ntt16> inv_phi_rom{};

        TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};

        explicit TwRom(int n) :
                N(n), w_rom(n >> 1), phi_rom(n), inv_w_rom(n >> 1), inv_phi_rom(n) {};

        static void initTwRom(TwRom &twRom, const int n) {
            twRom = TwRom(n);
        }
    };

    struct TwParam {
        std::vector<std::vector<Ntt16>> tw_factor{};

        TwParam() : tw_factor() {};

        explicit TwParam(int n) :
                tw_factor(n, std::vector<Ntt16>()) {};

        static void initTwParam(TwParam &twParam, const int n) {
            twParam = TwParam(n);
        }
    };

    extern TwParam NWC_TW;
    extern TwParam NWC_ITW;
    extern TwRom TW_ROM;


    void genTW_ROM(TwRom &tw_rom);

    void genNWCparam(TwParam &nwc_tw, const int n, const TwRom &tw_rom, const std::string &str);

    void applyNtt(Ntt16Polynomial &RES, const IntPolynomial &IN);

    void applyIntt(IntPolynomial &RES, const Ntt16Polynomial &IN);

    Ntt16 POW(Ntt16 BASE, Ntt16 EXP);

    Ntt16 modINV(Ntt16 in);

    Ntt16 modADD(Ntt16 a, Ntt16 b);

    Ntt16 modADDscale(Ntt16 a, Ntt16 b);

    Ntt16 modSUBscale(Ntt16 a, Ntt16 b);

    Ntt16 modSUB(Ntt16 a, Ntt16 b);

    Ntt16 redc(uint64_t T);

    Ntt16 modMULTR(Ntt16 a, Ntt16 b);

    Ntt16 modMULT(Ntt16 a, Ntt16 b);

    void initGlobalParamsNtt(int N);

    template<typename T, typename R>
    void applyNttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyNtt(out.a[row], in.a[row]);
        }
        applyNtt(out.b, in.b);
    }
}

#endif //YATFHE_NTT16_H
