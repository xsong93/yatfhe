//
// Created by Xintong Song on 2024/5/28.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/numeric_functions.h"

using namespace std;


TwParam16 NWC_TW16;
TwParam16 NWC_ITW16;
TwRom16 TW_ROM16;

Ntt16 POW16(Ntt16 BASE, Ntt16 EXP) {
    Ntt16 result = 1;
    while (EXP > 0) {
        if (EXP % 2 == 1) {
            result = uint32_t(result * BASE) % MOD16;
        }
        BASE = uint32_t(BASE * BASE) % MOD16;
        EXP = EXP / 2;
    }
    return result;
}

Ntt16 modINV16(Ntt16 in) {
    int32_t t = 0, newT = 1;
    int32_t r = MOD16, newR = in;
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
        throw std::invalid_argument("input is not invertible");
    }
    if (t < 0) {
        t += MOD16;
    }
    return static_cast<uint16_t>(t);
}

Ntt16 modADD16(Ntt16 a, Ntt16 b) {
    return ((MOD16 - a) > b) ? (a + b) : (a + b - MOD16);
}

Ntt16 modADDscale16(Ntt16 a, Ntt16 b) {
    Ntt16 temp = ((MOD16 - a) > b) ? (a + b) : (a + b - MOD16);

    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp >> 1) + ((MOD16 + 1) >> 1);
    }
    return temp;
}

Ntt16 modSUB16(Ntt16 a, Ntt16 b) {
    return (a >= b) ? (a - b) : (MOD16 - b + a);
}

Ntt16 modSUBscale16(Ntt16 a, Ntt16 b){
    Ntt16 temp = (a >= b) ? (a - b) : (MOD16 - b + a);
    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp = (temp >> 1) + ((MOD16 + 1) >> 1);
    }
    return temp;
}

Ntt16 modMULT16(Ntt16 a, Ntt16 b) {
    uint32_t result = static_cast<uint32_t>(a) * static_cast<uint32_t>(b);
    return static_cast<uint16_t>(result % MOD16);
}

void genTW_ROM16(TwRom16& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    Ntt16 w_q = Ntt16((MOD16 - 1) / (w_n << 1));
    Ntt16 phi_q = Ntt16((MOD16 - 1) / (phi_n << 1));
    Ntt16 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW16(PRIM_ROOT16, Ntt16(i * w_q));
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV16(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW16(PRIM_ROOT16, Ntt16(j*phi_q));
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV16(temp);
    }
}

void genNWCparam16(TwParam16& nwc_tw,const int n, const TwRom16& tw_rom, const std::string& str) {
    int lvl = calLogBase2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    Ntt16 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
    auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
    if (str == "NWC-DIT-NR-NNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = n >> (lvl - i);
            phi_probe = n >> (i + 1);
            scale = w_n >> i;
            for (int j = 0; j < tw_size; j++) {
                tw_temp = tw_rom.w_rom[j*scale];
                phi_temp = tw_rom.phi_rom[phi_probe];
                nwc_temp = modMULT16(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
            }
            bitRev(nwc_tw.tw_factor[i]);
        }
    }
    if (str == "NWC-DIF-RN-INNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = 1 << (lvl - i - 1);
            phi_probe = 1 << i;
            scale = 1 << i;
            for (int j = 0; j < tw_size; j++) {
                tw_temp = tw_rom.inv_w_rom[j*scale];
                phi_temp = tw_rom.inv_phi_rom[phi_probe];
                nwc_temp = modMULT16(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
            }
            bitRev(nwc_tw.tw_factor[i]);
        }
    }

}

void DIT_NR16(Ntt16Polynomial& RES, const Ntt16Polynomial& IN) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;

    const auto& tw = NWC_TW16.tw_factor;
    const auto lvl = calLogBase2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt16 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                temp_mult = modMULT16(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD16(res[j*block_size + k], temp_mult);
                temp_sub = modSUB16(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void applyNtt16(Ntt16Polynomial& RES, const IntPolynomial& IN) {
    auto N = IN.N;
    Ntt16Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt16(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt16(IN.coeffs[i] + MOD16);
        }
    }
    DIT_NR16(RES, format_input);
}

void DIF_RN16(Ntt16Polynomial& RES, const Ntt16Polynomial& IN) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = NWC_ITW16.tw_factor;
    const auto &N = IN.N;
    const auto lvl = calLogBase2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
//    Ntt16 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    Ntt16 temp_add, temp_sub, temp_mult;
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
                temp_add = modADDscale16(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_sub = modSUBscale16(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_mult = modMULT16(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void applyIntt16(IntPolynomial & RES, const Ntt16Polynomial& IN) {
    auto N = IN.N;
    Ntt16Polynomial res(N);
    DIF_RN16(res,IN);
    int16_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD16) {
            temp_ntt = int16_t(res.coeffs[i] - MOD16);
        } else {
            temp_ntt = int16_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT16_MASK);
        if (temp_poly >= POLY_MAX8) {
            RES.coeffs[i] = int16_t(temp_poly - (POLY_MAX8 << 1));
        } else {
            RES.coeffs[i] = int16_t(temp_poly);
        }
    }
}

void initGlobalParamsNtt16(int N) {
    auto depth = calLogBase2(N);
    TwParam16::initTwParam(NWC_TW16, depth);
    TwParam16::initTwParam(NWC_ITW16, depth);
    TwRom16::initTwRom(TW_ROM16, N);
    genTW_ROM16(TW_ROM16);
    genNWCparam16(NWC_TW16, N, TW_ROM16, STR_NTT);
    genNWCparam16(NWC_ITW16, N, TW_ROM16, STR_INTT);
}