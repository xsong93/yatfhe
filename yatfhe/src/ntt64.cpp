//
// Created by ic on 24-4-29.
//
#include "yatfhe/ntt64.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"
#include <iostream>
#include <gmp.h>
#include <string>
#include <cmath>

ntt64 POW(ntt64 BASE, ntt64 EXP, ntt64 MODU) {
    mpz_t base, exp, modu, res;
    ntt64 RES;
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

ntt64 modINV(ntt64 in){
    mpz_t a, inv, modu;
    mpz_inits(a, inv, modu, NULL);
    mpz_set_ui(a,in);
    mpz_set_ui(modu,MOD64);
    mpz_invert(inv,a,modu);
    ntt64 res = mpz_get_ui(inv);
    return res;
}

ntt64 modADD(ntt64 a, ntt64 b) {
    ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);
    return temp;
}

ntt64 modADDscale(ntt64 a, ntt64 b) {
    ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);

    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return ntt64(temp);
}
ntt64 modSUB(ntt64 a, ntt64 b) {
    ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    return temp;
}


ntt64 modSUBscale(ntt64 a, ntt64 b){
    ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return ntt64(temp);
}

ntt64 modMULT(ntt64 a, ntt64 b) {
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
    return ntt64(mpz_get_ui(TEMP));
}

void bit_rev(std::vector<ntt64>& x) {
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

void genTW_ROM(TW_ROM& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    ntt64 w_q = ntt64((MOD64 - 1)/(w_n<<1));
    ntt64 phi_q = ntt64((MOD64 - 1)/(phi_n<<1));
    ntt64 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW(PRIM_ROOT64, ntt64(i * w_q), MOD64);
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW(PRIM_ROOT64, ntt64(j*phi_q),MOD64);
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV(temp);
    }
}

void genNWCparam(TW_PARAM& nwc_tw,const int n, const TW_ROM& tw_rom, const std::string str) {
    int lvl = clog2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    ntt64 w_q = ntt64((MOD64 - 1)/(w_n<<1));
    ntt64 phi_q = ntt64((MOD64 - 1)/(phi_n<<1));
    ntt64 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
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
                nwc_temp = modMULT(tw_temp, phi_temp);
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
                nwc_temp = modMULT(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);
            }
            bit_rev(nwc_tw.tw_factor[i]);
        }
    }

}


void DIT_NR(Ntt64Polynomial& RES, const Ntt64Polynomial& IN, const TW_PARAM& ntt_param) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& N = IN.N;
    res = in;
//    std::vector<ntt64> in(N,0);

    const auto& tw = ntt_param.tw_factor;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    ntt64 temp_add, temp_sub, temp_mult, pos_a, pos_b;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (lvl-i);
        block_size = N >> i;
        gap = block_size >> 1;
        pos_a = 0; pos_b = 0;
        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                flag_a = res[j*block_size + k];
                flag_b = res[j*block_size + k + gap];
                flag_tw = tw[i][tw_index];
                pos_a = j*block_size + k;
                pos_b = j*block_size + k + gap;
                temp_mult = modMULT(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD(res[j*block_size + k], temp_mult);
                temp_sub = modSUB(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}
void doNTT(Ntt64Polynomial& RES, const IntPolynomial& IN, const TW_PARAM& ntt_param) {
    auto N = IN.N;
    Ntt64Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = ntt64(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = ntt64(IN.coeffs[i] + MOD64);
        }
    }
    DIT_NR(RES, format_input, ntt_param);
}


void DIF_RN(Ntt64Polynomial& RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = intt_param.tw_factor;
    const auto &N = IN.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    ntt64 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    ntt64 temp_add, temp_sub, temp_mult, pos_a, pos_b;
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

void doINTT(IntPolynomial & RES, const Ntt64Polynomial& IN, const TW_PARAM& intt_param) {
    auto N = IN.N;
    Ntt64Polynomial res(N);
    DIF_RN(res,IN,intt_param);
    int64_t temp_ntt = 0;
    int32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & mask);
        if (temp_poly >= poly_max) {
            RES.coeffs[i] = int32_t(temp_poly - poly_max);
        } else {
            RES.coeffs[i] = int32_t(temp_poly);
        }
    }
}

