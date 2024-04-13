//
// Created by ic on 24-4-10.
//
#include "yatfhe/myNtt.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"
#include <iostream>
#include <gmp.h>
#include <cmath>


void doNTT32 (NttPolynomial& res, const NttPolynomial& in, const Ntt32_TW& TW_param) {
    const auto N = in.N;
    auto lvl = clog2(in.N);
    auto gap = 0;
    auto tw_index = 0;
    auto interval = 0;
    auto& out = res.coeffs;
    const auto& tw = TW_param.tw_factor;
    Ntt32 tw_flag = 0;
    Ntt32 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < (N>>1); i++) {
        temp_add = modADD(in.coeffs[i], in.coeffs[i+(N>>1)]);
        temp_sub = modSUB(in.coeffs[i], in.coeffs[i+(N>>1)]);
//        tw_flag = tw[i];
        temp_mult = modMULT(temp_sub, tw[i]);
        out[i] = temp_add;
        out[i + (N>>1)] = temp_mult;
    }
    for (auto j = 1; j < lvl; j++) {
        gap = N>>(j+1);
        interval = 1<<j;
//        std::cout<<"gap = "<<gap<<"; interval = "<<interval<<std::endl;
        for (auto k = 0; k < interval; k++) {
            for (auto l = 0; l < gap ; l++) {
//                std::cout<<"k:l="<<k<<","<<l<<std::endl;
                tw_index = 1<<j;//
//                tw_flag = tw[l*tw_index];
                temp_add = modADD(out[k*gap*2 + l], out[k*gap*2 + gap + l]);
                temp_sub = modSUB(out[k*gap*2 + l], out[k*gap*2 + gap + l]);
                temp_mult = modMULT(temp_sub, tw[l*tw_index]);
                out[k*gap*2 + l] = temp_add;
                out[k*gap*2 + gap + l] = temp_mult;
            }
//            std::cout<<std::endl;
        }
    }
}


void doINTT32 (NttPolynomial& res, const NttPolynomial& in, const Ntt32_iTW& iTW_param) {
    const auto N = in.N;
    auto lvl = clog2(in.N);
    auto gap = 0;
    auto itw_index = 0;
    auto block = 0;
    auto& out = res.coeffs;
    const auto& itw = iTW_param.itw_factor;
    Ntt32 temp_add, temp_sub, temp_mult;
    for (auto i = 0; i < N; i = i+2) {
        Ntt32 debug1 = in.coeffs[i+1], debug2 = itw[0];
        temp_mult = modMULT(in.coeffs[i+1], itw[0]);
        temp_add = modADD(in.coeffs[i], temp_mult);
        temp_sub = modSUB(in.coeffs[i], temp_mult);

        out[i] = temp_add;
        out[i + 1] = temp_sub;
    }
    for (auto j = 1; j < lvl; j++) {
        gap = 1<<(j);
        block = N>>(j+1);
//        std::cout<<"gap = "<<gap<<"; block = "<<block<<std::endl;
        for (auto k = 0; k < block; k++) {
            for (auto l = 0; l < gap ; l++) {
//                std::cout<<"k:l="<<k<<","<<l<<std::endl;
                itw_index = N>>(j+1);
                temp_mult = modMULT(out[k*gap*2 + gap + l], itw[l*itw_index]);
                temp_add = modADD(out[k*gap*2 + l], temp_mult);
                temp_sub = modSUB(out[k*gap*2 + l], temp_mult);

                out[k*gap*2 + l] = temp_add;
                out[k*gap*2 + gap + l] = temp_sub;
            }
//            std::cout<<std::endl;
        }
    }
}
void print_myNtt(const Ntt32_TW& TW, const Ntt32_iTW& iTW) {
    std::cout<<"modulus is "<<MOD<<std::endl;
    std::cout<<"tw_factor: Q is "<<TW.Q<<std::endl;
    std::cout<<"itw_factor: Q is "<<iTW.Q<<std::endl;
    for (auto i = 0; i < TW.N; i++) {
        std::cout<<i<<":"<<TW.tw_factor[i]<<" ";
    }
    std::cout<<std::endl;
    for (auto i = 0; i < iTW.N; i++) {
        std::cout<<i<<":"<<iTW.itw_factor[i]<<" ";
    }
    std::cout<<std::endl;

}

void printNttPoly(const NttPolynomial& in) {
    for (auto i = 0; i < in.N; i++) {
        std::cout<<i<<":"<<in.coeffs[i]<<" ";
    }
    std::cout<<std::endl;
}

Ntt32 POW(Ntt32 BASE, Ntt32 EXP, Ntt32 MODU) {
    mpz_t base, exp, modu, res;
    uint32_t RES;
//    mpz
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

Ntt32 modINV(Ntt32 in){
    mpz_t a, inv, modu;
    mpz_inits(a, inv, modu, NULL);
    mpz_set_ui(a,in);
    mpz_set_ui(modu,MOD);
    mpz_invert(inv,a,modu);
    Ntt32 res = mpz_get_ui(inv);
    return res;
}

void genTW(Ntt32_TW& TW) {
    auto N = TW.N;
    for (auto i = 0; i < N ; i++) {
        TW.tw_factor[i] = POW(PRIM_ROOT,Ntt32(i*TW.Q), MOD);
    }
}
void geniTW(Ntt32_iTW& iTW, const Ntt32_TW& TW) {
    auto N = TW.N;
    for (auto i = 0; i < N; i++) {
        iTW.itw_factor[i] = modINV(TW.tw_factor[i]);
    }
}

//----------------------------------------------------------------------------

void genNTT32_PARAM(Ntt32_PARAM& ntt_param) {
    auto tw_n = ntt_param.tw_N;
    auto phi_n = ntt_param.phi_N;
    for (auto i = 0; i < tw_n; i++) {
        ntt_param.tw_factor[i] = POW(PRIM_ROOT, Ntt32(), MOD);
    }
}








//----------------------------------------------------------------------------



Ntt32 modADD(Ntt32 a, Ntt32 b) {
    int64_t temp;
    temp = int64_t(a) + int64_t(b);
    temp = temp >= MOD ? temp - MOD : temp;
    return Ntt32(temp);
}
Ntt32 modADDscale(Ntt32 a, Ntt32 b, bool isINTT) {
    int64_t temp;
    temp = int64_t(a) + int64_t(b);
    if (isINTT){
        temp = (temp>>1) >= MOD ? (temp>>1) - MOD : temp>>1;
    } else {
        temp = temp >= MOD ? temp - MOD : temp;
    }
    return Ntt32(temp);
}
Ntt32 modSUB(Ntt32 a, Ntt32 b) {
    int64_t temp = 0;
    temp = int64_t(a) - int64_t(b);
    temp = temp >= 0 ? temp : temp + MOD;
    return Ntt32(temp);
}
Ntt32 modSUBscale(Ntt32 a, Ntt32 b, bool isINTT){
    int64_t temp = 0;
    temp = (int64_t(a) - int64_t(b));
    if (isINTT) {
        temp = temp >= 0 ? temp/2 : temp/2 + MOD;
    } else {
        temp = temp >= 0 ? temp : temp + MOD;
    }
    return Ntt32(temp);
}

Ntt32 modMULT(Ntt32 a, Ntt32 b) {
    mpz_t A, B, TEMP, P;
    mpz_init(A);
    mpz_init(B);
    mpz_init(TEMP);
    mpz_init(P);
    mpz_set_ui(A,a);
    mpz_set_ui(B,b);
    mpz_set_ui(P,MOD);
    mpz_mul(TEMP, A, B);
    mpz_mod(TEMP, TEMP, P);
    return mpz_get_ui(TEMP);
}

int clog2(int N) {
    int res = 0;
    while (N >>= 1){
        res ++;
    }
    return res;
}



