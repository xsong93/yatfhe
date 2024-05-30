//
// Created by ic on 24-4-29.
//
#include <gmp.h>
#include <string>
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"

TwParam NWC_TW64;
TwParam NWC_ITW64;
TwRom TW_ROM64;

Ntt64 POW64(Ntt64 BASE, Ntt64 EXP, Ntt64 MODU) {
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

Ntt64 modINV64(Ntt64 in){
    mpz_t a, inv, modu;
    mpz_inits(a, inv, modu, NULL);
    mpz_set_ui(a,in);
    mpz_set_ui(modu,MOD64);
    mpz_invert(inv,a,modu);
    Ntt64 res = mpz_get_ui(inv);
    return res;
}

Ntt64 modADD64(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);
    return temp;
}

Ntt64 modADDscale64(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = ((MOD64 - a) > b) ? (a + b) : (a + b - MOD64);

    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return Ntt64(temp);
}
Ntt64 modSUB64(Ntt64 a, Ntt64 b) {
    Ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    return temp;
}


Ntt64 modSUBscale64(Ntt64 a, Ntt64 b){
    Ntt64 temp = 0;
    temp = (a >= b) ? (a - b) : (MOD64 - b + a);
    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp =  (temp>>1) + ((MOD64 + 1)>>1);
    }
    return Ntt64(temp);
}

Ntt64 modMULT64(Ntt64 a, Ntt64 b) {
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

// debug
Ntt64 modmul64(Ntt64 x, Ntt64 y) {
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

//https://1drv.ms/o/s!AvAq0B6vaN_Tj02p-JtAH0A9ICGw?e=ITjEAn
Ntt64 fastmm(Ntt64 x, Ntt64 y) {
    uint32_t x0 = (uint32_t)x;
    uint32_t x1 = (uint32_t)(x >> 32);
    uint32_t y0 = (uint32_t)y;
    uint32_t y1 = (uint32_t)(y >> 32);

    uint64_t x0y0 = (uint64_t)x0 * (uint64_t)y0;
    uint64_t x0y1 = (uint64_t)x0 * (uint64_t)y1;
    uint64_t x1y0 = (uint64_t)x1 * (uint64_t)y0;
    uint64_t x1y1 = (uint64_t)x1 * (uint64_t)y1;
    uint32_t x0y0_l = (uint32_t)x0y0;
    uint32_t x0y0_h = (uint32_t)(x0y0>>32);
    uint32_t x0y1_l = (uint32_t)x0y1;
    uint32_t x0y1_h = (uint32_t)(x0y1>>32);
    uint32_t x1y0_l = (uint32_t)x1y0;
    uint32_t x1y0_h = (uint32_t)(x1y0>>32);
    uint32_t x1y1_l = (uint32_t)x1y1;
    uint32_t x1y1_h = (uint32_t)(x1y1>>32);

    // z = x*y, z<127:0>
    // a = z<127:96>, b = z<95:64>, c = z<63:32>, d = z<31:0>
    uint32_t d = x0y0_l;
    uint64_t temp_c = (uint64_t)x0y0_h + (uint64_t)x0y1_l + (uint64_t)x1y0_l;
    uint32_t c = (uint32_t)temp_c;
    uint64_t c_of = (uint64_t)(temp_c >> 32);
    uint64_t temp_b = (uint64_t)x1y1_l + (uint64_t)x0y1_h + (uint64_t)x1y0_h + c_of;
    uint32_t b = (uint32_t)temp_b;
    uint32_t b_of = (uint32_t)(temp_b >> 32);
    uint32_t a = x1y1_h + b_of;
    uint64_t tmp_sum_bc = (uint64_t)c + (uint64_t)b;
    bool sumbc_of = (bool)(tmp_sum_bc>>32);
    tmp_sum_bc = (uint64_t) ((uint32_t)tmp_sum_bc)<<32;
    int64_t d_minus_ab = (int64_t)d - (int64_t)a - (int64_t)b;
    bool minus_flag = (d_minus_ab < 0);
//    uint64_t abs_val = abs(d_minus_ab);
    uint64_t res;
    if (sumbc_of) {
        if (minus_flag) {
//          res = (tmp_sum_bc >= abs_val)? tmp_sum_bc - abs_val + UINT64_MAX + 1 - MOD64 : UINT64_MAX + 1 + tmp_sum_bc - abs_val;
          res = (tmp_sum_bc >= (-d_minus_ab))? tmp_sum_bc + d_minus_ab + UINT64_MAX + 1 - MOD64 : UINT64_MAX + 1 + tmp_sum_bc + d_minus_ab;
          res = (res >= MOD64)? (res - MOD64):res;
        } else {
            res = tmp_sum_bc + d_minus_ab + UINT64_MAX + 1 - MOD64;
        }
    } else {
        if (minus_flag) {
            res = (tmp_sum_bc >= (-d_minus_ab))? tmp_sum_bc + d_minus_ab : MOD64 + tmp_sum_bc + d_minus_ab;
        } else {
            res = ((tmp_sum_bc + d_minus_ab) >= MOD64)? tmp_sum_bc + d_minus_ab - MOD64 : tmp_sum_bc + d_minus_ab;
        }
    }
    return res;

}

//    uint32_t tmp_l0 = (uint32_t)x0y1;
//    uint32_t tmp_l1 = (uint32_t)(x0y1>>32);
//    uint32_t tmp_h0 = (uint32_t)x1y0;
//    uint32_t tmp_h1 = (uint32_t)(x1y0>>32);
//    uint64_t tmp_add_l = tmp_l0 + tmp_l1;
//    uint64_t tmp_add_h = tmp_h0 + tmp_h1;
//    uint8_t tmp_overflow_l = (uint8_t)(tmp_add_l>>32);
//    uint8_t tmp_overflow_h = (uint8_t)((tmp_add_h + tmp_overflow_l)>>32);
//    uint64_t temp = (uint64_t)((uint32_t)tmp_add_l) + (uint64_t)((uint32_t)((tmp_add_h + tmp_overflow_l)>>32));

int clog2(int N) {
    int res = 0;
    while (N >>= 1){
        res ++;
    }
    return res;
}

void genTW_ROM64(TwRom& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    Ntt64 w_q = Ntt64((MOD64 - 1)/(w_n<<1));
    Ntt64 phi_q = Ntt64((MOD64 - 1)/(phi_n<<1));
    Ntt64 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW64(PRIM_ROOT64, Ntt64(i * w_q), MOD64);
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV64(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW64(PRIM_ROOT64, Ntt64(j*phi_q),MOD64);
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV64(temp);
    }
}

void genNWCparam64(TwParam& nwc_tw,const int n, const TwRom& tw_rom, const std::string& str) {
    int lvl = clog2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
//    Ntt64 w_q = Ntt64((MOD64 - 1)/(w_n<<1));
//    Ntt64 phi_q = Ntt64((MOD64 - 1)/(phi_n<<1));
    Ntt64 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
    auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
//    auto debug_tw = 0;
    if (str == "NWC-DIT-NR-NNT") {
        for (int i = 0; i < lvl; i++) {
            tw_size = n >> (lvl - i);
            phi_probe = n >> (i + 1);
            scale = w_n >> i;
            for (int j = 0; j < tw_size; j++) {
//                debug_tw = j*scale;
                tw_temp = tw_rom.w_rom[j*scale];
                phi_temp = tw_rom.phi_rom[phi_probe];
                nwc_temp = modMULT64(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);

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
//                debug_tw = j*scale;
                tw_temp = tw_rom.inv_w_rom[j*scale];
                phi_temp = tw_rom.inv_phi_rom[phi_probe];
                nwc_temp = modMULT64(tw_temp, phi_temp);
                nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);
            }
            bitRev(nwc_tw.tw_factor[i]);
        }
    }

}


void DIT_NR64(Ntt64Polynomial& RES, const Ntt64Polynomial& IN, const TwParam& ntt_param) {
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
                temp_mult = modMULT64(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD64(res[j*block_size + k], temp_mult);
                temp_sub = modSUB64(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void DIT_NR64(Ntt64Polynomial& RES, const Ntt64Polynomial& IN) {
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
                temp_mult = modMULT64(res[j*block_size + k + gap],tw[i][tw_index]);
                temp_add = modADD64(res[j*block_size + k], temp_mult);
                temp_sub = modSUB64(res[j*block_size + k], temp_mult);

                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_sub;
            }
        }
    }
}

void doNTT64(Ntt64Polynomial& RES, const IntPolynomial& IN, const TwParam& ntt_param) {
    auto N = IN.N;
    Ntt64Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt64(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt64(IN.coeffs[i] + MOD64);
        }
    }
    DIT_NR64(RES, format_input, ntt_param);
}

void doNTT64(Ntt64Polynomial& RES, const IntPolynomial& IN) {
    auto N = IN.N;
    Ntt64Polynomial format_input(N);
    for (int i = 0; i < N; i++) {
        if (IN.coeffs[i] >= 0){
            format_input.coeffs[i] = Ntt64(IN.coeffs[i]);
        } else {
            format_input.coeffs[i] = Ntt64(IN.coeffs[i] + MOD64);
        }
    }
    DIT_NR64(RES, format_input);
}


void DIF_RN64(Ntt64Polynomial& RES, const Ntt64Polynomial& IN, const TwParam& intt_param) {
    auto &res = RES.coeffs;
    const auto &in = IN.coeffs;
    const auto &tw = intt_param.tw_factor;
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
                temp_add = modADDscale64(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_sub = modSUBscale64(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_mult = modMULT64(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void DIF_RN64(Ntt64Polynomial& RES, const Ntt64Polynomial& IN) {
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
                temp_add = modADDscale64(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_sub = modSUBscale64(res[j * block_size + k], res[j * block_size + k + gap]);
                temp_mult = modMULT64(temp_sub, tw[i][tw_index]);
                res[j * block_size + k] = temp_add;
                res[j * block_size + k + gap] = temp_mult;
            }
        }
    }
}

void doINTT64(IntPolynomial & RES, const Ntt64Polynomial& IN, const TwParam& intt_param) {
    auto N = IN.N;
    Ntt64Polynomial res(N);
    DIF_RN64(res,IN,intt_param);
    int64_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT64_MASK);
        if (temp_poly >= POLY_MAX32) {
            RES.coeffs[i] = int32_t(temp_poly - (POLY_MAX32<<1));
        } else {
            RES.coeffs[i] = int32_t(temp_poly);
        }
    }
}

void doINTT64(IntPolynomial & RES, const Ntt64Polynomial& IN) {
    auto N = IN.N;
    Ntt64Polynomial res(N);
    DIF_RN64(res,IN);
    int64_t temp_ntt = 0;
    uint32_t temp_poly = 0;
    for (int i = 0; i < N; i++) {
        if (res.coeffs[i] >= HALF_MOD64) {
            temp_ntt = int64_t(res.coeffs[i] - MOD64);
        } else {
            temp_ntt = int64_t(res.coeffs[i]);
        }
        temp_poly = uint32_t(temp_ntt & NTT64_MASK);
        if (temp_poly >= POLY_MAX32) {
            RES.coeffs[i] = int32_t(temp_poly - (POLY_MAX32<<1));
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
    genTW_ROM64(TW_ROM64);
    genNWCparam64(NWC_TW64, N, TW_ROM64, STR_NTT);
    genNWCparam64(NWC_ITW64, N, TW_ROM64, STR_INTT);
}