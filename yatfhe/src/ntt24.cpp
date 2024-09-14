//
// Created by Xintong Song on 2024/9/13.
//

#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/numeric_functions.h"

using namespace std;

TwParam24 NWC_TW24;
TwParam24 NWC_ITW24;
TwRom24 TW_ROM24;

Ntt24 POW24(Ntt24 BASE, Ntt24 EXP) {
    Ntt24 result = 1;
    uint64_t t = BASE;
    while (EXP > 0) {
        if (EXP % 2 == 1) {
            result = static_cast<uint64_t>(result * t) % MOD24;
        }
        t = static_cast<uint64_t>(t * t) % MOD24;
        EXP = EXP / 2;
    }
    return result;
}

Ntt24 modINV24(Ntt24 in) {
    int32_t t = 0;
    int32_t newT = 1;
    int32_t r = MOD24;
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
        t += MOD24;
    }
    return static_cast<Ntt24>(t);
}

Ntt24 modADD24(Ntt24 a, Ntt24 b) {
    return ((MOD24 - a) > b) ? (a + b) : (a + b - MOD24);
}

Ntt24 modADDscale24(Ntt24 a, Ntt24 b) {
    Ntt24 temp = ((MOD24 - a) > b) ? (a + b) : (a + b - MOD24);

    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp >> 1) + ((MOD24 + 1) >> 1);
    }
    return temp;
}

Ntt24 modSUB24(Ntt24 a, Ntt24 b) {
    return (a >= b) ? (a - b) : (MOD24 - b + a);
}

Ntt24 modSUBscale24(Ntt24 a, Ntt24 b){
    Ntt24 temp = (a >= b) ? (a - b) : (MOD24 - b + a);
    if (temp % 2 == 0) {
        temp = temp >> 1;
    } else {
        temp = (temp >> 1) + ((MOD24 + 1) >> 1);
    }
    return temp;
}

Ntt24 modMULT24(Ntt24 a, Ntt24 b) {
    uint64_t a1 = a;
    uint64_t b1 = b;
    return static_cast<Ntt24>(a1 * b1 % MOD24);
}

void genTW_ROM24(TwRom24& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    auto w_q = static_cast<Ntt24>((MOD24 - 1) / (w_n << 1));
    auto phi_q = static_cast<Ntt24>((MOD24 - 1) / (phi_n << 1));
    Ntt24 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW24(PRIM_ROOT24, static_cast<Ntt24>(i * w_q));
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV24(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW24(PRIM_ROOT24, static_cast<Ntt24>(j*phi_q));
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV24(temp);
    }
}

void genNWCparam24(TwParam24& nwc_tw,const int n, const TwRom24& tw_rom, const std::string& str) {
    int lvl = calLogBase2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    Ntt24 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
    auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
    if (str == "NWC-DIT-NR-NNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = n >> (lvl - i);
            phi_probe = n >> (i + 1);
            scale = w_n >> i;
            for (int j = 0; j < tw_size; j++) {
                tw_temp = tw_rom.w_rom[j*scale];
                phi_temp = tw_rom.phi_rom[phi_probe];
                nwc_temp = modMULT24(tw_temp, phi_temp);
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
                nwc_temp = modMULT24(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
            }
            bitRev(nwc_tw.tw_factor[i]);
        }
    }

}

void DIT_NR24(Ntt24Polynomial& RES, const Ntt24Polynomial& IN) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;

    const auto& tw = NWC_TW24.tw_factor;
    const auto lvl = calLogBase2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt24 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                temp_mult = modMULT24(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD24(res[j*block_size + k], temp_mult);
                temp_sub = modSUB24(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void applyNtt24(Ntt24Polynomial& RES, const Int8Polynomial& IN) {
    auto N = IN.N;
    Ntt24Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = static_cast<Ntt24>(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = static_cast<Ntt24>(IN.coeffs[i] + MOD24);
        }
    }
    DIT_NR24(RES, format_input);
}

void DIF_RN24(Ntt24Polynomial& RES, const Ntt24Polynomial& IN) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = NWC_ITW24.tw_factor;
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
                temp_add = modADDscale24(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_sub = modSUBscale24(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_mult = modMULT24(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void applyIntt24(Int8Polynomial & RES, const Ntt24Polynomial& IN) {
    auto N = IN.N;
    Ntt24Polynomial res(N);
    DIF_RN24(res,IN);
    int32_t temp_ntt;
    uint8_t temp_poly;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD24) {
            temp_ntt = static_cast<int32_t>(res.coeffs[i] - MOD24);
        } else {
            temp_ntt = static_cast<int32_t>(res.coeffs[i]);
        }
        temp_poly = static_cast<uint8_t>(temp_ntt & NTT24_MASK);
        if (temp_poly >= POLY_MAX8) {
            RES.coeffs[i] = static_cast<int8_t>(temp_poly - (POLY_MAX8 << 1));
        } else {
            RES.coeffs[i] = static_cast<int8_t>(temp_poly);
        }
    }
}

void initGlobalParamsNtt24(int N) {
    auto depth = calLogBase2(N);
    TwParam24::initTwParam(NWC_TW24, depth);
    TwParam24::initTwParam(NWC_ITW24, depth);
    TwRom24::initTwRom(TW_ROM24, N);
    genTW_ROM24(TW_ROM24);
    genNWCparam24(NWC_TW24, N, TW_ROM24, STR_NTT);
    genNWCparam24(NWC_ITW24, N, TW_ROM24, STR_INTT);
}

void modularAccumulate24(std::vector<Ntt24>& coeffsB, const std::vector<Ntt24>& coeffsA, const std::vector<Ntt24>& coeffsS) {
    const auto N = coeffsB.size();
    for (auto j = 0; j < N; j++) {
        auto tmp = modMULT24(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modADD24(coeffsB[j], tmp);
    }
}

void calModularInnerProductNtt24(Ntt24Polynomial& b, const Ntt24Polynomial& a, const Ntt24Polynomial& s) {
    modularAccumulate24(b.coeffs, a.coeffs, s.coeffs);
}