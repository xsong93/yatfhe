//
// Created by ic on 24-4-10.
//
#include "yatfhe/myNtt.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/polynomial.h"
#include <iostream>
#include <gmp.h>
#include <string>
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


Ntt32 modADD(Ntt32 a, Ntt32 b) {
    int64_t temp;
    temp = int64_t(a) + int64_t(b);
    temp = temp >= MOD ? temp - MOD : temp;
    return Ntt32(temp);
}
Ntt32 modADDscale(Ntt32 a, Ntt32 b) {
    int64_t temp = 0;
    temp = (int64_t(a) + int64_t(b));
    temp = (temp) >= MOD ? (temp - MOD) : temp;
    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp = ( temp + MOD + 1) >> 1;
    }
    return Ntt32(temp);
}
Ntt32 modSUB(Ntt32 a, Ntt32 b) {
    int64_t temp = 0;
    temp = int64_t(a) - int64_t(b);
    temp = (temp < 0)? temp + MOD : temp;
    return Ntt32(temp);
}

Ntt32 modSUBscale(Ntt32 a, Ntt32 b){
    int64_t temp = 0;
    temp = (int64_t(a) - int64_t(b)); // debug: where minus is negative, add p first or do the div first?
    temp = (temp < 0)? temp + MOD : temp;
    if (temp%2 == 0) {
        temp = temp >> 1;
    } else {
        temp = (temp + MOD + 1) >> 1;
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


//----------------------------------------------------------------------------



void genTW_ROM(TW_ROM& tw_rom) {
    auto w_n = (tw_rom.N) >> 1;
    auto phi_n = tw_rom.N;
    Ntt32 w_q = Ntt32((MOD - 1)/(w_n<<1));
    Ntt32 phi_q = Ntt32((MOD - 1)/(phi_n<<1));
    Ntt32 temp = 0;
    for (int i = 0; i < w_n; i++) {
        temp = POW(PRIM_ROOT, Ntt32(i * w_q), MOD);
        tw_rom.w_rom[i] = temp;
        tw_rom.inv_w_rom[i] = modINV(temp);
    }
    for (int j = 0; j < phi_n; j++) {
        temp = POW(PRIM_ROOT, Ntt32(j*phi_q),MOD);
        tw_rom.phi_rom[j] = temp;
        tw_rom.inv_phi_rom[j] = modINV(temp);
    }

}
void genNTT32_PARAM(Ntt32_PARAM& ntt_param) {
    auto tw_n = ntt_param.tw_N;
    auto phi_n = ntt_param.phi_N;
    Ntt32 tw_q = Ntt32((MOD-1)/ntt_param.N);
    Ntt32 phi_q = Ntt32((MOD - 1)/(ntt_param.N << 1));
    for (auto i = 0; i < tw_n; i++) {
        ntt_param.tw_factor[i] = POW(PRIM_ROOT, Ntt32(i*tw_q), MOD);
    }
    for (auto j = 0; j < phi_n; j++) {
        ntt_param.phi_factor[j] = POW(PRIM_ROOT, Ntt32(j*phi_q),MOD);
    }
}
void genROM(ROM& rom){
    auto w_n = rom.N;
    auto phi_n = rom.N * 2;
    Ntt32 w_q = Ntt32((MOD - 1)/(w_n<<1));
    Ntt32 phi_q = Ntt32((MOD - 1)/(phi_n<<1));
    for (auto i = 0; i < w_n; i++) {
        rom.ntt_rom.tw_factor[i] = POW(PRIM_ROOT, Ntt32(i*w_q), MOD);
        rom.intt_rom.inv_tw_factor[i] = modINV(rom.ntt_rom.tw_factor[i]);
    }
    bit_rev(rom.intt_rom.inv_tw_factor);
    for (auto j = 0; j < phi_n; j++) {
        rom.ntt_rom.phi_factor[j] = POW(PRIM_ROOT, j*phi_q, MOD);
        rom.intt_rom.inv_phi_factor[j] = modINV(rom.ntt_rom.phi_factor[j]);
    }
}

void genNWCparam(TW_PARAM& nwc_tw,const int n, const TW_ROM& tw_rom, const std::string str) {
    int lvl = clog2(n);
    auto w_n = n >> 1;
    auto phi_n = n;
    Ntt32 w_q = Ntt32((MOD - 1)/(w_n<<1));
    Ntt32 phi_q = Ntt32((MOD - 1)/(phi_n<<1));
    Ntt32 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
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


void genDIF_ROM(DIF_ROM& rom) {
    const auto lvl = rom.l;
    auto w_n = 1 << (lvl-1);
    auto phi_n = 1 << lvl;
    auto N = 1 << lvl;
    Ntt32 w_q = Ntt32((MOD - 1)/(w_n<<1));
    Ntt32 phi_q = Ntt32((MOD - 1)/(phi_n<<1));
    Ntt32 temp = 0, temp_inv = 0, scale = 0;
    auto phi_size = 0, tw_size = 0, probe = 0;
    for (int i = 0; i < lvl; i++) {
        scale = 1 << i;
        phi_size = phi_n >> i;
        tw_size = w_n >> i;
        for (int j = 0; j < tw_size; j++) {
            probe = j*scale;
            temp = POW(PRIM_ROOT, Ntt32(j*w_q*scale), MOD);
            temp_inv = modINV(temp);
            rom.ntt_tw.tw_factor[i].push_back(temp);
            rom.intt_tw.tw_factor[i].push_back(temp_inv);
        }
        bit_rev(rom.intt_tw.tw_factor[i]);
//        for (int k = 0; k < phi_size; k++) {
//            temp = POW(PRIM_ROOT, Ntt32(k*phi_q*scale), MOD);
//            temp_inv = modINV(temp);
//            rom.phi_tw.tw_factor[i].push_back(temp);
//            rom.iphi_tw.tw_factor[i].push_back(temp_inv);
//        }
//        bit_rev(rom.iphi_tw.tw_factor[i]);
    }
    temp = 0;


}


void pre_process(NttPolynomial& in, const Ntt32_PARAM& para) {
    auto N = para.phi_N;
    auto& res = in.coeffs;
    auto& coeff = para.phi_factor;
    for (auto i = 0; i < N; i++) {
        res[i] = modMULT(res[i], coeff[i]);
    }
}


void DIT_NR(NttPolynomial& RES, const NttPolynomial& IN, const TW_PARAM& ntt_param) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& tw = ntt_param.tw_factor;
    const auto& N = IN.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt32 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    Ntt32 temp_add, temp_sub, temp_mult, pos_a, pos_b;
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

void DIF_NR(NttPolynomial& RES, const NttPolynomial& IN, const Ntt32_PARAM& ntt_param) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& tw = ntt_param.tw_factor;
    const auto& N = IN.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    Ntt32 temp_add, temp_sub, temp_mult;
    res = in;
    for (auto i = 0; i < lvl; i++) {
         block = 1 << i;
         block_size = N >> i;
         gap = N >> (i+1);
        for (auto j = 0; j < block; j++) {
            for (auto k = 0; k < gap; k++) {
                temp_add = modADD(res[j*block_size + k], res[j*block_size + k + gap]);
                temp_sub = modSUB(res[j*block_size + k], res[j*block_size + k + gap]);
                temp_mult = modMULT(temp_sub,tw[k*block]);
                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_mult;
            }
        }

    }
}



void DIF_RN(NttPolynomial& RES, const NttPolynomial& IN, const TW_PARAM& intt_param) {
    auto& res = RES.coeffs;
    const auto& in = IN.coeffs;
    const auto& tw = intt_param.tw_factor;
    const auto& N = IN.N;
    const auto lvl = clog2(N);
    auto gap = 0;
    auto block = 0;
    auto block_size = 0;
    auto tw_index = 0;
    Ntt32 flag_a = 0, flag_b = 0, flag_tw = 0;
    res = in;
    Ntt32 temp_add, temp_sub, temp_mult, pos_a, pos_b;
    for (auto i = 0; i < lvl; i++) {
        block = N >> (i + 1);
        block_size = 1 << (i + 1);
        gap = 1 << i;
        pos_a = 0; pos_b = 0;

        for (auto j = 0; j < block; j++) { //debug:tw_index overflow
            tw_index = j;
            for (auto k = 0; k < gap; k++) {
                flag_a = res[j*block_size + k];
                flag_b = res[j*block_size + k + gap];
                flag_tw = tw[i][tw_index];
                pos_a = j*block_size + k;
                pos_b = j*block_size + k + gap;
                temp_add = modADDscale(res[j*block_size + k], res[j*block_size + k + gap]);
                temp_sub = modSUBscale(res[j*block_size + k], res[j*block_size + k + gap]);
                temp_mult = modMULT(temp_sub,tw[i][tw_index]);
                res[j*block_size + k] = temp_add;
                res[j*block_size + k + gap] = temp_mult;

            }
        }
    }
}







//----------------------------------------------------------------------------

void printROM(ROM& rom) {
    auto w_n = rom.N;
    auto phi_n = rom.N * 2;
    std::cout<<"ntt rom:"<<std::endl;
    std::cout<<"ntt tw:"<<std::endl;
    for (auto i = 0; i < w_n; i++) {
        std::cout<<i<<":"<<rom.ntt_rom.tw_factor[i]<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"ntt phi:"<<std::endl;
    for (auto j = 0; j < phi_n; j++) {
        std::cout << j << ":" << rom.ntt_rom.phi_factor[j] << " ";
    }
    std::cout<<std::endl;
    std::cout<<"intt rom:"<<std::endl;
    std::cout<<"intt inv_tw:"<<std::endl;
    for (auto k = 0; k < w_n; k++) {
        std::cout << k << ":" << rom.intt_rom.inv_tw_factor[k] << " ";
    }
    std::cout<<std::endl;
    std::cout<<"intt inv_phi:"<<std::endl;
    for (auto l = 0; l < phi_n; l++) {
        std::cout << l << ":" << rom.intt_rom.inv_phi_factor[l] << " ";
    }
    std::cout<<std::endl;
}


void bit_rev(std::vector<Ntt32>& x) {
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
//void bitreverse (std::vector<Ntt32>& x, int N) {
//    int j = 0;
//    int b = 0;
//
//    for (int i = 1; i < N; i++) {
//        b = N >> 1;  // Initialize b to half of N
//        while (j >= b) {
//            j -= b;  // Perform bit-reversal
//            b >>= 1;
//        }
//        j += b;  // Move to the next position
//
//        // Swap elements if the bit-reversed index is greater than the current index
//        if (j > i) {
//            NttType temp = x[j];
//            x[j] = x[i];
//            x[i] = temp;
//        }
//    }
//}
//
