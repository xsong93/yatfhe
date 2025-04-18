//
// Created by Xintong Song on 2024/9/13.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/numeric_functions.h"

using namespace std;
using namespace NttNative;

namespace NttNative24 {

    TwParam NWC_TW;
    TwParam NWC_ITW;
    TwRom TW_ROM;

    Ntt24 POW(Ntt24 BASE, Ntt24 EXP) {
        Ntt24 result = 1;
        uint64_t t = BASE;
        while (EXP > 0) {
            if (EXP % 2 == 1) {
                result = static_cast<uint64_t>(result * t) % MOD;
            }
            t = static_cast<uint64_t>(t * t) % MOD;
            EXP = EXP / 2;
        }
        return result;
    }

    Ntt24 modINV(Ntt24 in) {
        int32_t t = 0;
        int32_t newT = 1;
        int32_t r = MOD;
        auto newR = static_cast<int32_t>(in);
        while (newR != 0) {
            int32_t quotient = r / newR;
            int32_t tempT = newT;
            newT = t - quotient * newT;
            t = tempT;

            int32_t tempR = newR;
            newR = r - quotient * newR;
            r = tempR;
        }
        if (t < 0) {
            t += MOD;
        }
        return static_cast<Ntt24>(t);
    }

    Ntt24 modADD(Ntt24 a, Ntt24 b) {
        return ((MOD - a) > b) ? (a + b) : (a + b - MOD);
    }

    Ntt24 modADDscale(Ntt24 a, Ntt24 b) {
        Ntt24 temp = ((MOD - a) > b) ? (a + b) : (a + b - MOD);

        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt24 modSUB(Ntt24 a, Ntt24 b) {
        return (a >= b) ? (a - b) : (MOD - b + a);
    }

    Ntt24 modSUBscale(Ntt24 a, Ntt24 b) {
        Ntt24 temp = (a >= b) ? (a - b) : (MOD - b + a);
        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt24 modMULT(Ntt24 a, Ntt24 b) {
        uint64_t a1 = a;
        uint64_t b1 = b;
        return static_cast<Ntt24>(a1 * b1 % MOD);
    }

    void genTW_ROM(TwRom &tw_rom) {
        auto w_n = (tw_rom.N) >> 1;
        auto phi_n = tw_rom.N;
        auto w_q = static_cast<Ntt24>((MOD - 1) / (w_n << 1));
        auto phi_q = static_cast<Ntt24>((MOD - 1) / (phi_n << 1));
        Ntt24 temp = 0;
        for (int i = 0; i < w_n; i++) {
            temp = POW(PRIM_ROOT, static_cast<Ntt24>(i * w_q));
            tw_rom.w_rom[i] = temp;
            tw_rom.inv_w_rom[i] = modINV(temp);
        }
        for (int j = 0; j < phi_n; j++) {
            temp = POW(PRIM_ROOT, static_cast<Ntt24>(j * phi_q));
            tw_rom.phi_rom[j] = temp;
            tw_rom.inv_phi_rom[j] = modINV(temp);
        }
    }

    void genNWCparam(TwParam &nwc_tw, int n, const TwRom &tw_rom, const std::string &str) {
        int lvl = calLogBase2(n);
        auto w_n = n >> 1;
        auto phi_n = n;
        Ntt24 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
        auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
        if (str == "NWC-DIT-NR-NTT") {
            for (int i = 0; i < lvl; i++) {
                tw_size = n >> (lvl - i);
                phi_probe = n >> (i + 1);
                scale = w_n >> i;
                for (int j = 0; j < tw_size; j++) {
                    tw_temp = tw_rom.w_rom[j * scale];
                    phi_temp = tw_rom.phi_rom[phi_probe];
                    nwc_temp = modMULT(tw_temp, phi_temp);
                    nwc_tw.tw_factor[i].push_back(nwc_temp);
                }
                bitRev(nwc_tw.tw_factor[i]);
            }
        }
        if (str == "NWC-DIF-RN-INTT") {
            for (int i = 0; i < lvl; i++) {
                tw_size = 1 << (lvl - i - 1);
                phi_probe = 1 << i;
                scale = 1 << i;
                for (int j = 0; j < tw_size; j++) {
                    tw_temp = tw_rom.inv_w_rom[j * scale];
                    phi_temp = tw_rom.inv_phi_rom[phi_probe];
                    nwc_temp = modMULT(tw_temp, phi_temp);
                    nwc_tw.tw_factor[i].push_back(nwc_temp);
                }
                bitRev(nwc_tw.tw_factor[i]);
            }
        }

    }

    void DIT_NR24(Ntt24Polynomial &RES, const Ntt24Polynomial &IN) {
        auto &res = RES.coeffs;
        const auto &in = IN.coeffs;
        const auto &N = IN.N;
        res = in;

        const auto &tw = NWC_TW.tw_factor;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
        Ntt24 temp_add, temp_sub, temp_mult;
        for (auto i = 0; i < lvl; i++) {
            block = N >> (lvl - i);
            block_size = N >> i;
            gap = block_size >> 1;
            for (auto j = 0; j < block; j++) { //debug:tw_index overflow
                tw_index = j;
                for (auto k = 0; k < gap; k++) {
                    temp_mult = modMULT(res[j * block_size + k + gap], tw[i][tw_index]);
                    temp_add = modADD(res[j * block_size + k], temp_mult);
                    temp_sub = modSUB(res[j * block_size + k], temp_mult);

                    res[j * block_size + k] = temp_add;
                    res[j * block_size + k + gap] = temp_sub;
                }
            }
        }
    }

    void applyNtt(Ntt24Polynomial &RES, const Int8Polynomial &IN) {
        auto N = IN.N;
        Ntt24Polynomial format_input{N};
        for (int i = 0; i < N; i++) {
            auto &valIn = IN.coeffs[i];
            if (valIn >= 0) {
                format_input.coeffs[i] = static_cast<Ntt24>(valIn);
            } else {
                format_input.coeffs[i] = static_cast<Ntt24>(valIn + MOD);
            }
        }
        DIT_NR24(RES, format_input);
    }

    void DIF_RN24(Ntt24Polynomial &RES, const Ntt24Polynomial &IN) {
        auto &res = RES.coeffs;
        const auto &in = IN.coeffs;
        const auto &tw = NWC_ITW.tw_factor;
        const auto &N = IN.N;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
//    Ntt24 flag_a = 0, flag_b = 0, flag_tw = 0;
        res = in;
        Ntt24 temp_add, temp_sub, temp_mult;
//    int pos_a, pos_b;
        for (auto i = 0; i < lvl; i++) {
            block = N >> (i + 1);
            block_size = 1 << (i + 1);
            gap = 1 << i;
//        pos_a = 0;
//        pos_b = 0;
            for (auto j = 0; j < block; j++) { //debug:tw_index overflow
                tw_index = j;
                for (auto k = 0; k < gap; k++) {
//                flag_a = res[j * block_size + k];
//                flag_b = res[j * block_size + k + gap];
//                flag_tw = tw[i][tw_index];
//                pos_a = j * block_size + k;
//                pos_b = j * block_size + k + gap;
                    temp_add = modADDscale(res[j * block_size + k], res[j * block_size + k + gap]);
                    temp_sub = modSUBscale(res[j * block_size + k], res[j * block_size + k + gap]);
                    temp_mult = modMULT(temp_sub, tw[i][tw_index]);
                    res[j * block_size + k] = temp_add;
                    res[j * block_size + k + gap] = temp_mult;
                }
            }
        }
    }

    void applyIntt(Int8Polynomial &out, const Ntt24Polynomial &in, int q) {
        auto N = in.N;
        Ntt24Polynomial res(N);
        DIF_RN24(res, in);
        int32_t temp_ntt;
        int32_t temp_poly;
        int qHalf = q / 2;
        int loHalf = -qHalf;
        int hiHalf = (q % 2 == 0) ? (qHalf - 1) : qHalf;
        for (int i = 0; i < N; i++) {
            if (res.coeffs[i] >= HALF_MOD) {
                temp_ntt = static_cast<int32_t>(res.coeffs[i] - MOD);
            } else {
                temp_ntt = static_cast<int32_t>(res.coeffs[i]);
            }
//        temp_poly = static_cast<uint8_t>(temp_ntt & NTT_MASK);
//        if (temp_poly >= POLY_MAX8) {
//            RES.coeffs[i] = static_cast<int8_t>(temp_poly - (POLY_MAX8 << 1));
//        } else {
//            RES.coeffs[i] = static_cast<int8_t>(temp_poly);
//        }
            temp_poly = temp_ntt % q;
            if (temp_poly < loHalf) {
                out.coeffs[i] = static_cast<int8_t>(temp_poly + q);
            } else if (temp_poly > hiHalf) {
                out.coeffs[i] = static_cast<int8_t>(temp_poly - q);
            } else {
                out.coeffs[i] = static_cast<int8_t>(temp_poly);
            }
        }
    }

    void initGlobalParamsNtt(int N) {
        auto depth = calLogBase2(N);
        TwParam::initTwParam(NWC_TW, depth);
        TwParam::initTwParam(NWC_ITW, depth);
        TwRom::initTwRom(TW_ROM, N);
        genTW_ROM(TW_ROM);
        genNWCparam(NWC_TW, N, TW_ROM, STR_NTT);
        genNWCparam(NWC_ITW, N, TW_ROM, STR_INTT);
    }

    void modularAccumulate24(std::vector<Ntt24> &coeffsOut, const std::vector<Ntt24> &coeffs1,
                             const std::vector<Ntt24> &coeffs2) {
        const auto N = coeffsOut.size();
        for (auto j = 0; j < N; j++) {
            auto tmp = modMULT(coeffs1[j], coeffs2[j]);
            coeffsOut[j] = modADD(coeffsOut[j], tmp);
        }
    }

    void calModularInnerProductNtt(Ntt24Polynomial &out, const Ntt24Polynomial &in1, const Ntt24Polynomial &in2) {
        modularAccumulate24(out.coeffs, in1.coeffs, in2.coeffs);
    }
}