//
// Created by Xintong Song on 20325/3/4.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt32.h"
#include "yatfhe/numeric.h"

using namespace std;
using namespace NttNative;

namespace NttNative32 {

    TwParam NWC_TW;
    TwParam NWC_ITW;
    TwRom TW_ROM;

    Ntt32 POW(Ntt32 BASE, Ntt32 EXP) {
        Ntt32 result = 1;
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

    Ntt32 modINV(Ntt32 in) {
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
        return static_cast<Ntt32>(t);
    }

    Ntt32 modADD(Ntt32 a, Ntt32 b) {
        return ((MOD - a) > b) ? (a + b) : (a + b - MOD);
    }

    Ntt32 modADDscale(Ntt32 a, Ntt32 b) {
        Ntt32 temp = ((MOD - a) > b) ? (a + b) : (a + b - MOD);

        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt32 modSUB(Ntt32 a, Ntt32 b) {
        return (a >= b) ? (a - b) : (MOD - b + a);
    }

    Ntt32 modSUBscale(Ntt32 a, Ntt32 b) {
        Ntt32 temp = (a >= b) ? (a - b) : (MOD - b + a);
        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt32 modMULT(Ntt32 a, Ntt32 b) {
        uint64_t a1 = a;
        uint64_t b1 = b;
        return static_cast<Ntt32>(a1 * b1 % MOD);
    }

    void genTW_ROM(TwRom &tw_rom) {
        auto w_n = (tw_rom.N) >> 1;
        auto phi_n = tw_rom.N;
        auto w_q = static_cast<Ntt32>((MOD - 1) / (w_n << 1));
        auto phi_q = static_cast<Ntt32>((MOD - 1) / (phi_n << 1));
        Ntt32 temp = 0;
        for (int i = 0; i < w_n; i++) {
            temp = POW(PRIM_ROOT, static_cast<Ntt32>(i * w_q));
            tw_rom.w_rom[i] = temp;
            tw_rom.inv_w_rom[i] = modINV(temp);
        }
        for (int j = 0; j < phi_n; j++) {
            temp = POW(PRIM_ROOT, static_cast<Ntt32>(j * phi_q));
            tw_rom.phi_rom[j] = temp;
            tw_rom.inv_phi_rom[j] = modINV(temp);
        }
    }

    void genNWCparam(TwParam &nwc_tw, int n, const TwRom &tw_rom, const std::string &str) {
        int lvl = calLogBase2(n);
        auto w_n = n >> 1;
        auto phi_n = n;
        Ntt32 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
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

    void DIT_NR32(Ntt32Polynomial &RES, const Ntt32Polynomial &IN) {
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
        Ntt32 temp_add, temp_sub, temp_mult;
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

    void applyNtt(Ntt32Polynomial &RES, const TorusPolynomial &IN) {
        auto N = IN.N;
        Ntt32Polynomial format_input{N};
        for (int i = 0; i < N; i++) {
            auto &valIn = IN.coeffs[i];
            if (valIn >= 0) {
                format_input.coeffs[i] = static_cast<Ntt32>(valIn);
            } else {
                format_input.coeffs[i] = static_cast<Ntt32>(valIn + MOD);
            }
        }
        DIT_NR32(RES, format_input);
    }

    void DIF_RN32(Ntt32Polynomial &RES, const Ntt32Polynomial &IN) {
        auto &res = RES.coeffs;
        const auto &in = IN.coeffs;
        const auto &tw = NWC_ITW.tw_factor;
        const auto &N = IN.N;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
//    Ntt32 flag_a = 0, flag_b = 0, flag_tw = 0;
        res = in;
        Ntt32 temp_add, temp_sub, temp_mult;
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

    void applyIntt(TorusPolynomial &out, const Ntt32Polynomial &in) {
        auto N = in.N;
        Ntt32Polynomial res(N);
        DIF_RN32(res, in);
        for (int i = 0; i < N; i++) {
            if (res.coeffs[i] >= HALF_MOD) {
                out.coeffs[i] = static_cast<int32_t>(res.coeffs[i] - MOD);
            } else {
                out.coeffs[i] = static_cast<int32_t>(res.coeffs[i]);
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

    void modularAccumulate32(std::vector<Ntt32> &coeffsOut, const std::vector<Ntt32> &coeffs1,
                             const std::vector<Ntt32> &coeffs2) {
        const auto N = coeffsOut.size();
        for (auto j = 0; j < N; j++) {
            auto tmp = modMULT(coeffs1[j], coeffs2[j]);
            coeffsOut[j] = modADD(coeffsOut[j], tmp);
        }
    }

    void calModularInnerProductNtt(Ntt32Polynomial &out, const Ntt32Polynomial &in1, const Ntt32Polynomial &in2) {
        modularAccumulate32(out.coeffs, in1.coeffs, in2.coeffs);
    }
}