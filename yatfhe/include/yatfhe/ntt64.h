//
// Created by ic on 24-4-29.
//

#ifndef HLS_YATFHE_NTT64_H
#define HLS_YATFHE_NTT64_H

#include <vector>
#include <string>
#include <gmp.h>
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"

namespace NttNative64 {

    constexpr Ntt64 MOD = 0xffffffff00000001;
    constexpr Ntt64 HALF_MOD = (MOD + 1) >> 1;
    constexpr uint32_t NTT_MASK = 0xffffffff;
    constexpr Ntt64 PRIM_ROOT = 7;
    constexpr uint32_t POLY_MAX = 1 << 31;

    struct TwRom {
        int N{};
        std::vector<Ntt64> w_rom{};
        std::vector<Ntt64> inv_w_rom{};
        std::vector<Ntt64> phi_rom{};
        std::vector<Ntt64> inv_phi_rom{};

        TwRom() : N(), w_rom(), phi_rom(), inv_w_rom(), inv_phi_rom() {};

        explicit TwRom(int n) :
                N(n), w_rom(n >> 1), phi_rom(n), inv_w_rom(n >> 1), inv_phi_rom(n) {};

        static void initTwRom(TwRom &twRom, const int n) {
            twRom = TwRom(n);
        }
    };

    struct TwParam {
        std::vector<std::vector<Ntt64>> tw_factor{};

        TwParam() : tw_factor() {};

        explicit TwParam(int n) :
                tw_factor(n, std::vector<Ntt64>()) {};

        static void initTwParam(TwParam &twParam, const int n) {
            twParam = TwParam(n);
        }
    };


    extern TwParam NWC_TW;
    extern TwParam NWC_ITW;
    extern TwRom TW_ROM;

//----------------------------------------------------------------------------------


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

    void genTW_ROM(TwRom &tw_rom);

    void genNWCparam(TwParam &nwc_tw, const int n, const TwRom &tw_rom, const std::string &str);

    Ntt64 POW(Ntt64 base, Ntt64 exp, Ntt64 mod);

    Ntt64 modINV(Ntt64 in);

    Ntt64 modAdd(Ntt64 x, Ntt64 y);

    Ntt64 modADDscale(Ntt64 a, Ntt64 b);

    Ntt64 modSUBscale(Ntt64 a, Ntt64 b);

    Ntt64 modSub(Ntt64 x, Ntt64 y);

    Ntt64 modMULT(Ntt64 a, Ntt64 b);

    Ntt64 fastmm(Ntt64 a, Ntt64 b);

    Ntt64 fastmm_opt(Ntt64 a, Ntt64 b);

    Ntt64 modmul64(Ntt64 x, Ntt64 y);

    void applyNtt(NttPolynomial &out, const IntPolynomial &in);

    void applyIntt(IntPolynomial &out, const NttPolynomial &in);

    void
    modularMult(std::vector<NttType> &output, const std::vector<NttType> &coeffsA, const std::vector<NttType> &coeffsB);

    void modularAccumulate(vector<NttType> &res, const vector<NttType> &in1, const vector<NttType> &in2);

    void calModularInnerProductNtt(NttPolynomial &res, const vector<NttPolynomial> &in1,
                                   const vector<NttPolynomial> &in2);

    void
    calModularInnerProductNtt(NttPolynomial &res, const NttPolynomial &in1, const NttPolynomial &in2);

    void initGlobalParamsNtt(int N);
}

#endif //HLS_YATFHE_NTT64_H
