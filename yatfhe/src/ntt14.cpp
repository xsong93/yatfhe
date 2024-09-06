//
// Created by Xintong Song on 2024/5/28.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/numeric_functions.h"

using namespace std;

TwParam14 NWC_TW14;
TwParam14 NWC_ITW14;
TwRom14 TW_ROM14;

Ntt14 POW14(Ntt14 BASE, Ntt14 EXP) {
    Ntt14 result = 1;
    while (EXP > 0) {
        if (EXP % 2 == 1) {
            result = uint32_t(result * BASE) % MOD14;
        }
        BASE = uint32_t(BASE * BASE) % MOD14;
        EXP = EXP / 2;
    }
    return result;
}

Ntt14 modINV14(Ntt14 in) {
    int32_t t = 0;
    int32_t newT = 1;
    int32_t r = MOD14;
    int32_t newR = in;
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
        t += MOD14;
    }
    return static_cast<Ntt14>(t);
}

Ntt14 modADD14(Ntt14 a, Ntt14 b) {
    return ((MOD14 - a) > b) ? (a + b) : (a + b - MOD14);
}

Ntt14 modADDscale14(Ntt14 a, Ntt14 b) {
    Ntt14 temp = ((MOD14 - a) > b) ? (a + b) : (a + b - MOD14);

    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp >> 1) + ((MOD14 + 1) >> 1);
    }
    return temp;
}

Ntt14 modSUB14(Ntt14 a, Ntt14 b) {
    return (a >= b) ? (a - b) : (MOD14 - b + a);
}

Ntt14 modSUBscale14(Ntt14 a, Ntt14 b){
    Ntt14 temp = (a >= b) ? (a - b) : (MOD14 - b + a);
    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp = (temp >> 1) + ((MOD14 + 1) >> 1);
    }
    return temp;
}

Ntt14 modMULT14(Ntt14 a, Ntt14 b) {
    uint32_t result = a * b;
    return static_cast<Ntt14>(result % MOD14);
}

void genTW_ROM14(TwRom14& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    Ntt14 w_q = Ntt14((MOD14 - 1) / (w_n << 1));
    Ntt14 phi_q = Ntt14((MOD14 - 1) / (phi_n << 1));
    Ntt14 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW14(PRIM_ROOT14, Ntt14(i * w_q));
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV14(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW14(PRIM_ROOT14, Ntt14(j*phi_q));
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV14(temp);
    }
}

void genNWCparam14(TwParam14& nwc_tw,const int n, const TwRom14& tw_rom, const std::string& str) {
    int lvl = calLogBase2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    Ntt14 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
    auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
    if (str == "NWC-DIT-NR-NNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = n >> (lvl - i);
            phi_probe = n >> (i + 1);
            scale = w_n >> i;
            for (int j = 0; j < tw_size; j++) {
                tw_temp = tw_rom.w_rom[j*scale];
                phi_temp = tw_rom.phi_rom[phi_probe];
                nwc_temp = modMULT14(tw_temp, phi_temp);
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
                nwc_temp = modMULT14(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
            }
            bitRev(nwc_tw.tw_factor[i]);
        }
    }

}

void DIT_NR14(Ntt14Polynomial& RES, const Ntt14Polynomial& IN) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;

    const auto& tw = NWC_TW14.tw_factor;
    const auto lvl = calLogBase2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt14 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                temp_mult = modMULT14(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD14(res[j*block_size + k], temp_mult);
                temp_sub = modSUB14(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void applyNtt14(Ntt14Polynomial& RES, const IntPolynomial& IN) {
    auto N = IN.N;
    Ntt14Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt14(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt14(IN.coeffs[i] + MOD14);
        }
    }
    DIT_NR14(RES, format_input);
}

void applyNtt14Poly8(Ntt14Polynomial& RES, const Int8Polynomial& IN) {
    auto N = IN.N;
    Ntt14Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt14(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt14(IN.coeffs[i] + MOD14);
        }
    }
    DIT_NR14(RES, format_input);
}

void DIF_RN14(Ntt14Polynomial& RES, const Ntt14Polynomial& IN) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = NWC_ITW14.tw_factor;
    const auto &N = IN.N;
    const auto lvl = calLogBase2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt14 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    Ntt14 temp_add, temp_sub, temp_mult;
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
                temp_add = modADDscale14(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_sub = modSUBscale14(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_mult = modMULT14(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void applyIntt14(IntPolynomial & RES, const Ntt14Polynomial& IN) {
    auto N = IN.N;
    Ntt14Polynomial res(N);
    DIF_RN14(res,IN);
    int32_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD14) {
            temp_ntt = int32_t(res.coeffs[i] - MOD14);
        } else {
            temp_ntt = int32_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT14_MASK);
        if (temp_poly >= POLY_MAX8) {
            RES.coeffs[i] = int32_t(temp_poly - (POLY_MAX8 << 1));
        } else {
            RES.coeffs[i] = int32_t(temp_poly);
        }
    }
}

void initGlobalParamsNtt14(int N) {
    auto depth = calLogBase2(N);
    TwParam14::initTwParam(NWC_TW14, depth);
    TwParam14::initTwParam(NWC_ITW14, depth);
    TwRom14::initTwRom(TW_ROM14, N);
    genTW_ROM14(TW_ROM14);
    genNWCparam14(NWC_TW14, N, TW_ROM14, STR_NTT);
    genNWCparam14(NWC_ITW14, N, TW_ROM14, STR_INTT);
}

void modularAccumulate14(std::vector<Ntt14>& coeffsB, const std::vector<Ntt14>& coeffsA, const std::vector<Ntt14>& coeffsS) {
    const auto N = coeffsB.size();
    for (auto j = 0; j < N; j++) {
        auto tmp = modMULT14(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modADD14(coeffsB[j], tmp);
    }
}

void calModularInnerProductNtt14(Ntt14Polynomial& b, const Ntt14Polynomial& a, const Ntt14Polynomial& s) {
    modularAccumulate14(b.coeffs, a.coeffs, s.coeffs);
}