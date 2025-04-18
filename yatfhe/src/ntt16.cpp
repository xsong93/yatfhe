//
// Created by Xintong Song on 2024/5/28.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/numeric_functions.h"
#include <gmp.h>


using namespace std;
using namespace NttNative;

namespace NttNative16 {

    TwParam NWC_TW;
    TwParam NWC_ITW;
    TwRom TW_ROM;

    Ntt16 POW(Ntt16 BASE, Ntt16 EXP) {
        Ntt16 result = 1;
        while (EXP > 0) {
            if (EXP % 2 == 1) {
                result = uint32_t(result * BASE) % MOD;
            }
            BASE = uint32_t(BASE * BASE) % MOD;
            EXP = EXP / 2;
        }
        return result;
    }

    Ntt16 modinv16(Ntt16 in) {
        mpz_t a, inv, modu;
        mpz_init(a);
        mpz_init(inv);
        mpz_init(modu);
        mpz_set_ui(a, in);
        mpz_set_ui(modu, MOD);
        mpz_invert(inv, a, modu);
        Ntt16 res = mpz_get_ui(inv);
        return res;
    }

    Ntt16 modINV(Ntt16 in) {
        int32_t t = 0, newT = 1;
        int32_t r = MOD, newR = in;
        while (newR != 0) {
            int32_t quotient = r / newR;
            int32_t tempT = newT;
            newT = t - quotient * newT;
            t = tempT;

            int32_t tempR = newR;
            newR = r - quotient * newR;
            r = tempR;
        }
        if (r > 1) {
//        throw std::invalid_argument("input is not invertible");
        }
        if (t < 0) {
            t += MOD;
        }
        return static_cast<uint16_t>(t);
    }

    Ntt16 modADD(Ntt16 a, Ntt16 b) {
        return ((MOD - a) > b) ? (a + b) : (a + b - MOD);
    }

    Ntt16 modADDscale(Ntt16 a, Ntt16 b) {
        Ntt16 temp = ((MOD - a) > b) ? (a + b) : (a + b - MOD);

        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt16 modSUB(Ntt16 a, Ntt16 b) {
        return (a >= b) ? (a - b) : (MOD - b + a);
    }

    Ntt16 modSUBscale(Ntt16 a, Ntt16 b) {
        Ntt16 temp = (a >= b) ? (a - b) : (MOD - b + a);
        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return temp;
    }

    Ntt16 redc(uint64_t T) {
        T = (uint32_t) T + (T >> 32);
        return (T & 0xFFFF) + MOD - (T >> 16);
    }

    Ntt16 modMULTR(Ntt16 a, Ntt16 b) {
        return redc((uint64_t) a * b);
    }

    Ntt16 modMULT(Ntt16 a, Ntt16 b) {
        uint64_t result = a * b;
        return static_cast<uint32_t>(result % MOD);
    }

    void genTW_ROM(TwRom &tw_rom) {
        auto w_n = (tw_rom.N) >> 1;
        auto phi_n = tw_rom.N;
        Ntt16 w_q = Ntt16((MOD - 1) / (w_n << 1));
        Ntt16 phi_q = Ntt16((MOD - 1) / (phi_n << 1));
        Ntt16 temp = 0;
        for (int i = 0; i < w_n; i++) {
            temp = POW(PRIM_ROOT, Ntt16(i * w_q));
            tw_rom.w_rom[i] = temp;
            tw_rom.inv_w_rom[i] = modINV(temp);
        }
        for (int j = 0; j < phi_n; j++) {
            temp = POW(PRIM_ROOT, Ntt16(j * phi_q));
            tw_rom.phi_rom[j] = temp;
            tw_rom.inv_phi_rom[j] = modINV(temp);
        }
    }

    void genNWCparam(TwParam &nwc_tw, const int n, const TwRom &tw_rom, const std::string &str) {
        int lvl = calLogBase2(n);
        auto w_n = n >> 1;
        auto phi_n = n;
        Ntt16 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
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

    void DIT_NR16(Ntt16Polynomial &RES, const Ntt16Polynomial &IN) {
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
        Ntt16 temp_add, temp_sub, temp_mult;
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

    void applyNtt(Ntt16Polynomial &RES, const IntPolynomial &IN) {
        auto N = IN.N;
        Ntt16Polynomial format_input(N);
        for (int i = 0; i < N; i++) {
            if (IN.coeffs[i] >= 0) {
                format_input.coeffs[i] = Ntt16(IN.coeffs[i]);
            } else {
                format_input.coeffs[i] = Ntt16(IN.coeffs[i] + MOD);
            }
        }
        DIT_NR16(RES, format_input);
    }

    void DIF_RN16(Ntt16Polynomial &RES, const Ntt16Polynomial &IN) {
        auto &res = RES.coeffs;
        const auto &in = IN.coeffs;
        const auto &tw = NWC_ITW.tw_factor;
        const auto &N = IN.N;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
        Ntt16 flag_a = 0, flag_b = 0, flag_tw = 0;
        res = in;
        Ntt16 temp_add, temp_sub, temp_mult;
        int pos_a, pos_b;
        for (auto i = 0; i < lvl; i++) {
            block = N >> (i + 1);
            block_size = 1 << (i + 1);
            gap = 1 << i;
            pos_a = 0;
            pos_b = 0;
            for (auto j = 0; j < block; j++) { //debug:tw_index overflow
                tw_index = j;
                for (auto k = 0; k < gap; k++) {
                    flag_a = res[j * block_size + k];
                    flag_b = res[j * block_size + k + gap];
                    flag_tw = tw[i][tw_index];
                    pos_a = j * block_size + k;
                    pos_b = j * block_size + k + gap;
                    temp_add = modADDscale(res[j * block_size + k], res[j * block_size + k + gap]);
                    temp_sub = modSUBscale(res[j * block_size + k], res[j * block_size + k + gap]);
                    temp_mult = modMULT(temp_sub, tw[i][tw_index]);
                    res[j * block_size + k] = temp_add;
                    res[j * block_size + k + gap] = temp_mult;
                }
            }
        }
    }

    void applyIntt(IntPolynomial &RES, const Ntt16Polynomial &IN) {
        auto N = IN.N;
        Ntt16Polynomial res(N);
        DIF_RN16(res, IN);
        int32_t temp_ntt = 0;
        uint32_t temp_poly = 0;
        for (int i = 0; i < N; i++) {
            if (res.coeffs[i] >= HALF_MOD) {
                temp_ntt = int32_t(res.coeffs[i] - MOD);
            } else {
                temp_ntt = int32_t(res.coeffs[i]);
            }
            temp_poly = uint32_t(temp_ntt & NTT_MASK);
            if (temp_poly >= POLY_MAX8) {
                RES.coeffs[i] = int32_t(temp_poly - (POLY_MAX8 << 1));
            } else {
                RES.coeffs[i] = int32_t(temp_poly);
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
}