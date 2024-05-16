//
// Created by ic on 24-4-29.
//
#include <iostream>
#include <gmp.h>
#include <string>
#include <cmath>
#include "yatfhe/ntt64.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/yatfhe_parameters.h"

TwParam NWC_TW64;
TwParam NWC_ITW64;
TwRom TW_ROM64;

Ntt64 POW(Ntt64 BASE, Ntt64 EXP, Ntt64 MODU) {
    mpz_t base, exp, modu, res;
    Ntt64 RES;
    mpz_init(base);
    mpz_init(exp);
    mpz_init(modu);
    mpz_init(res);
    mpz_set_ui(base, BASE);
    mpz_set_ui(exp, EXP);
    mpz_set_ui(modu, MODU);
    mpz_powm(res, base, exp, modu);
    RES = mpz_get_ui(res);
    return RES;
}

Ntt64 modINV(Ntt64 in){
    mpz_t a, inv, modu;
    mpz_inits(a, inv, modu, NULL);
    mpz_set_ui(a,in);
    mpz_set_ui(modu,MOD64);
    mpz_invert(inv,a,modu);
    Ntt64 res = mpz_get_ui(inv);
    return res;
}

Ntt64 modADD(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);
    return temp;
}

Ntt64 modADDscale(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);

    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return Ntt64(temp);
}
Ntt64 modSUB(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    return temp;
}


Ntt64 modSUBscale(Ntt64 a, Ntt64 b){
    Ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return Ntt64(temp);
}

Ntt64 modMULT(Ntt64 a, Ntt64 b) {
    mpz_t A, B, TEMP, P;
    mpz_init(A);
    mpz_init(B);
    mpz_init(TEMP);
    mpz_init(P);
    mpz_set_ui(A,a);
    mpz_set_ui(B,b);
    mpz_set_ui(P,MOD64);
    mpz_mul(TEMP, A, B);
    mpz_mod(TEMP, TEMP, P);
    return Ntt64(mpz_get_ui(TEMP));
}

Ntt64 modmul(Ntt64 x, Ntt64 y) {
    auto x0 = (uint32_t)x;
    auto x1 = (uint32_t)(x >> 32);
    auto y0 = (uint32_t)y;
    auto y1 = (uint32_t)(y >> 32);

    // Perform 64-bit multiplication
    uint64_t x0y0 = (uint64_t)x0 * (uint64_t)y0;
    uint64_t x0y1 = (uint64_t)x0 * (uint64_t)y1;
    uint64_t x1y0 = (uint64_t)x1 * (uint64_t)y0;
    uint64_t x1y1 = (uint64_t)x1 * (uint64_t)y1;

    // Compute partial products and handle carry
    auto d = (uint32_t)x0y0;
    uint64_t pp1 = (x0y0 >> 32) + (uint32_t)(x1y0) + (uint32_t)(x0y1);
    auto c = (uint32_t)pp1;
    uint64_t pp2 = (x1y0 >> 32) + (x0y1 >> 32) + (uint32_t)(x1y1);
    uint64_t pp3 = (pp1 >> 32) + (uint32_t)(pp2);

    // Handle overflow and underflow
    uint32_t a = (pp2 >> 32) + (x1y1 >> 32);
    uint64_t bpc = (uint32_t)pp3 + (uint64_t)c;
    bpc = ((bpc + (bpc >> 32)) << 32) - (bpc >> 32);
    uint64_t minus = ((uint64_t)a + ((uint64_t)(uint32_t)pp3));
    uint64_t plus = bpc + (uint64_t)d;

    // Return the result modulo MODULUS
    if (plus >= minus) {
        return (plus - minus);
    }
    return MOD64 - minus + plus;
}

void bit_rev(std::vector<Ntt64>& x) {
    int j = 0;
    int b = 0;
    int N = int(x.size());
    for (int i = 1; i < N; i++) {
        b = N >> 1;  // Initialize b to half of N
        while (j >= b) {
            j -= b;  // Perform bit-reversal
            b >>= 1;
        }
        j += b;  // Move to the next position

        // Swap elements if the bit-reversed index is greater than the current index
        if (j > i) {
            NttType temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

int clog2(int N) {
    int res = 0;
    while (N >>= 1){
        res ++;
    }
    return res;
}

void genTW_ROM(TwRom& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    Ntt64 w_q = Ntt64((MOD64 - 1)/(w_n<<1));
    Ntt64 phi_q = Ntt64((MOD64 - 1)/(phi_n<<1));
    Ntt64 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW(PRIM_ROOT64, Ntt64(i * w_q), MOD64);
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW(PRIM_ROOT64, Ntt64(j*phi_q),MOD64);
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV(temp);
    }
}

void genNWCparam(TwParam& nwc_tw,const int n, const TwRom& tw_rom, const std::string str) {
    int lvl = clog2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    Ntt64 w_q = Ntt64((MOD64 - 1)/(w_n<<1));
    Ntt64 phi_q = Ntt64((MOD64 - 1)/(phi_n<<1));
    Ntt64 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
    auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
    auto debug_tw = 0;
    if (str == "NWC-DIT-NR-NNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = n >> (lvl - i);
            phi_probe = n >> (i + 1);
            scale = w_n >> i;
            for (int j = 0; j < tw_size; j++) {
                debug_tw = j*scale;
                tw_temp = tw_rom.w_rom[j*scale];
                phi_temp = tw_rom.phi_rom[phi_probe];
                nwc_temp = modmul(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);

            }
            bit_rev(nwc_tw.tw_factor[i]);
        }
    }
    if (str == "NWC-DIF-RN-INNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = 1 << (lvl - i - 1);
            phi_probe = 1 << i;
            scale = 1 << i;
            for (int j = 0; j < tw_size; j++) {
                debug_tw = j*scale;
                tw_temp = tw_rom.inv_w_rom[j*scale];
                phi_temp = tw_rom.inv_phi_rom[phi_probe];
                nwc_temp = modmul(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);
            }
            bit_rev(nwc_tw.tw_factor[i]);
        }
    }

}


void DIT_NR(Ntt64Polynomial& RES, const Ntt64Polynomial& IN, const TwParam& ntt_param) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;
//    std::vector<Ntt64> in(N,0);

    const auto& tw = ntt_param.tw_factor;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
//    Ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    Ntt64 temp_add, temp_sub, temp_mult;
//    int pos_a, pos_b;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
//        pos_a = 0; pos_b = 0;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
//                flag_a = res[j*block_size + k];
//                flag_b = res[j*block_size + k + gap];
//                flag_tw = tw[i][tw_index];
//                pos_a = j*block_size + k;
//                pos_b = j*block_size + k + gap;
                temp_mult = modmul(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD(res[j*block_size + k], temp_mult);
                temp_sub = modSUB(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void DIT_NR(Ntt64Polynomial& RES, const Ntt64Polynomial& IN) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;
//    std::vector<Ntt64> in(N,0);

    const auto& tw = NWC_TW64.tw_factor;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
//    Ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    Ntt64 temp_add, temp_sub, temp_mult;
//    int pos_a, pos_b;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
//        pos_a = 0; pos_b = 0;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
//                flag_a = res[j*block_size + k];
//                flag_b = res[j*block_size + k + gap];
//                flag_tw = tw[i][tw_index];
//                pos_a = j*block_size + k;
//                pos_b = j*block_size + k + gap;
                temp_mult = modmul(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD(res[j*block_size + k], temp_mult);
                temp_sub = modSUB(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void doNTT(Ntt64Polynomial& RES, const IntPolynomial& IN, const TwParam& ntt_param) {
    auto N = IN.N;
    Ntt64Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt64(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt64(IN.coeffs[i] + MOD64);
        }
    }
    DIT_NR(RES, format_input, ntt_param);
}

void doNTT(Ntt64Polynomial& RES, const IntPolynomial& IN) {
    auto N = IN.N;
    Ntt64Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt64(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt64(IN.coeffs[i] + MOD64);
        }
    }
    DIT_NR(RES, format_input, NWC_TW64);
}


void DIF_RN(Ntt64Polynomial& RES, const Ntt64Polynomial& IN) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = NWC_ITW64.tw_factor;
    const auto &N = IN.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
//    Ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    Ntt64 temp_add, temp_sub, temp_mult;
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
                temp_mult = modmul(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void doINTT(IntPolynomial & RES, const Ntt64Polynomial& IN) {
    auto N = IN.N;
    Ntt64Polynomial res(N);
    DIF_RN(res,IN);
    int64_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT64_MASK);
        if (temp_poly >= POLY_MAX) {
            RES.coeffs[i] = int32_t(temp_poly - (POLY_MAX<<1));
        } else {
            RES.coeffs[i] = int32_t(temp_poly);
        }
    }
}

void initGlobalParamsNtt64(int N) {
    auto depth = clog2(N);
    TwParam::initTwParam(NWC_TW64, depth);
    TwParam::initTwParam(NWC_ITW64, depth);
    TwRom::initTwRom(TW_ROM64, N);
    genTW_ROM(TW_ROM64);
    genNWCparam(NWC_TW64, N, TW_ROM64, STR_NTT);
    genNWCparam(NWC_ITW64, N, TW_ROM64, STR_INTT);
}
