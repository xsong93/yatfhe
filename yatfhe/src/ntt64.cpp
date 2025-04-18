//
// Created by ic on 24-4-29.
//
#include <gmp.h>
#include <string>
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"

using namespace NttNative;

namespace NttNative64 {

    TwParam NWC_TW;
    TwParam NWC_ITW;
    TwRom TW_ROM;

    void DITNRLaPoly(LagrangePolynomial &out, const LagrangePolynomial &in) {
        auto &output = out.coeffs;
        const auto &input = in.coeffs;
        const auto &N = in.N;
        output = input;

        const auto &tw = NWC_TW.tw_factor;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
        Ntt64 temp_add, temp_sub, temp_mult;
        for (auto i = 0; i < lvl; i++) {
            block = N >> (lvl - i);
            block_size = N >> i;
            gap = block_size >> 1;
            for (auto j = 0; j < block; j++) {
                tw_index = j;
                for (auto k = 0; k < gap; k++) {
                    temp_mult = fastmm_opt(output[j * block_size + k + gap], tw[i][tw_index]);
                    temp_add = modAdd(output[j * block_size + k], temp_mult);
                    temp_sub = modSub(output[j * block_size + k], temp_mult);

                    output[j * block_size + k] = temp_add;
                    output[j * block_size + k + gap] = temp_sub;
                }
            }
        }
    }

    void DIFRNLaPoly(LagrangePolynomial &out, const LagrangePolynomial &in) {
        auto &output = out.coeffs;
        const auto &input = in.coeffs;
        const auto &tw = NWC_ITW.tw_factor;
        const auto &N = in.N;
        const auto lvl = calLogBase2(N);
        auto gap = 0;
        auto block = 0;
        auto block_size = 0;
        auto tw_index = 0;
        output = input;
        Ntt64 temp_add, temp_sub, temp_mult;
//    Ntt64 tw_flag = 0;
        for (auto i = 0; i < lvl; i++) {
            block = N >> (i + 1);
            block_size = 1 << (i + 1);
            gap = 1 << i;
            for (auto j = 0; j < block; j++) { //debug:tw_index overflow
                tw_index = j;
                for (auto k = 0; k < gap; k++) {
//                tw_flag = tw[i][tw_index];
                    temp_add = modADDscale(output[j * block_size + k], output[j * block_size + k + gap]);
                    temp_sub = modSUBscale(output[j * block_size + k], output[j * block_size + k + gap]);
                    temp_mult = fastmm_opt(temp_sub, tw[i][tw_index]);
                    output[j * block_size + k] = temp_add;
                    output[j * block_size + k + gap] = temp_mult;
                }
            }
        }
    }

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

    Ntt64 modINV(Ntt64 in) {
        mpz_t a, inv, modu;
        mpz_inits(a, inv, modu, NULL);
        mpz_set_ui(a, in);
        mpz_set_ui(modu, MOD);
        mpz_invert(inv, a, modu);
        Ntt64 res = mpz_get_ui(inv);
        return res;
    }

    Ntt64 modAdd(Ntt64 x, Ntt64 y) {
        return ((MOD - x) > y) ? (x + y) : (x + y - MOD);
    }

    Ntt64 modADDscale(Ntt64 a, Ntt64 b) {
        Ntt64 temp = 0;
        temp = ((MOD - a) > b) ? (a + b) : (a + b - MOD);

        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return Ntt64(temp);
    }

    Ntt64 modSub(Ntt64 x, Ntt64 y) {
        return (x >= y) ? (x - y) : (MOD - y + x);
    }

    Ntt64 modSUBscale(Ntt64 a, Ntt64 b) {
        Ntt64 temp = 0;
        temp = (a >= b) ? (a - b) : (MOD - b + a);
        if (temp % 2 == 0) {
            temp = temp >> 1;
        } else {
            temp = (temp >> 1) + ((MOD + 1) >> 1);
        }
        return Ntt64(temp);
    }

    Ntt64 modMULT(Ntt64 a, Ntt64 b) {
        mpz_t A, B, TEMP, P;
        mpz_init(A);
        mpz_init(B);
        mpz_init(TEMP);
        mpz_init(P);
        mpz_set_ui(A, a);
        mpz_set_ui(B, b);
        mpz_set_ui(P, MOD);
        mpz_mul(TEMP, A, B);
        mpz_mod(TEMP, TEMP, P);
        return Ntt64(mpz_get_ui(TEMP));
    }

// debug
    Ntt64 modmul64(Ntt64 x, Ntt64 y) {
        auto x0 = (uint32_t) x;
        auto x1 = (uint32_t) (x >> 32);
        auto y0 = (uint32_t) y;
        auto y1 = (uint32_t) (y >> 32);

        // Perform 64-bit multiplication
        uint64_t x0y0 = (uint64_t) x0 * (uint64_t) y0;
        uint64_t x0y1 = (uint64_t) x0 * (uint64_t) y1;
        uint64_t x1y0 = (uint64_t) x1 * (uint64_t) y0;
        uint64_t x1y1 = (uint64_t) x1 * (uint64_t) y1;

        // Compute partial products and handle carry
        auto d = (uint32_t) x0y0;
        uint64_t pp1 = (x0y0 >> 32) + (uint32_t) (x1y0) + (uint32_t) (x0y1);
        auto c = (uint32_t) pp1;
        uint64_t pp2 = (x1y0 >> 32) + (x0y1 >> 32) + (uint32_t) (x1y1);
        uint64_t pp3 = (pp1 >> 32) + (uint32_t) (pp2);

        // Handle overflow and underflow
        uint32_t a = (pp2 >> 32) + (x1y1 >> 32);
        uint64_t bpc = (uint32_t) pp3 + (uint64_t) c;
        bpc = ((bpc + (bpc >> 32)) << 32) - (bpc >> 32);
        uint64_t minus = ((uint64_t) a + ((uint64_t) (uint32_t) pp3));
        uint64_t plus = bpc + (uint64_t) d;

        // Return the result modulo MODULUS
        if (plus >= minus) {
            return (plus - minus);
        }
        return MOD - minus + plus;
    }

//https://1drv.ms/o/s!AvAq0B6vaN_Tj02p-JtAH0A9ICGw?e=ITjEAn
    Ntt64 fastmm(Ntt64 x, Ntt64 y) {
        uint32_t x0 = (uint32_t) x;
        uint32_t x1 = (uint32_t) (x >> 32);
        uint32_t y0 = (uint32_t) y;
        uint32_t y1 = (uint32_t) (y >> 32);

        uint64_t x0y0 = (uint64_t) x0 * (uint64_t) y0;
        uint64_t x0y1 = (uint64_t) x0 * (uint64_t) y1;
        uint64_t x1y0 = (uint64_t) x1 * (uint64_t) y0;
        uint64_t x1y1 = (uint64_t) x1 * (uint64_t) y1;
        uint32_t x0y0_l = (uint32_t) x0y0;
        uint32_t x0y0_h = (uint32_t) (x0y0 >> 32);
        uint32_t x0y1_l = (uint32_t) x0y1;
        uint32_t x0y1_h = (uint32_t) (x0y1 >> 32);
        uint32_t x1y0_l = (uint32_t) x1y0;
        uint32_t x1y0_h = (uint32_t) (x1y0 >> 32);
        uint32_t x1y1_l = (uint32_t) x1y1;
        uint32_t x1y1_h = (uint32_t) (x1y1 >> 32);

        // z = x*y, z<127:0>
        // a = z<127:96>, b = z<95:64>, c = z<63:32>, d = z<31:0>
        uint32_t d = x0y0_l;
        uint64_t temp_c = (uint64_t) x0y0_h + (uint64_t) x0y1_l + (uint64_t) x1y0_l;
        uint32_t c = (uint32_t) temp_c;
        uint64_t c_of = (uint64_t) (temp_c >> 32);
        uint64_t temp_b = (uint64_t) x1y1_l + (uint64_t) x0y1_h + (uint64_t) x1y0_h + c_of;
        uint32_t b = (uint32_t) temp_b;
        uint32_t b_of = (uint32_t) (temp_b >> 32);
        uint32_t a = x1y1_h + b_of;
        uint64_t tmp_sum_bc = (uint64_t) c + (uint64_t) b;
        bool sumbc_of = (bool) (tmp_sum_bc >> 32);
        tmp_sum_bc = (uint64_t) ((uint32_t) tmp_sum_bc) << 32;
        int64_t d_minus_ab = (int64_t) d - (int64_t) a - (int64_t) b;
        bool minus_flag = (d_minus_ab < 0);
//    uint64_t abs_val = abs(d_minus_ab);
        uint64_t res;
        if (sumbc_of) {
            if (minus_flag) {
//          res = (tmp_sum_bc >= abs_val)? tmp_sum_bc - abs_val + UINT64_MAX + 1 - MOD : UINT64_MAX + 1 + tmp_sum_bc - abs_val;
                res = (tmp_sum_bc >= (-d_minus_ab)) ? tmp_sum_bc + d_minus_ab + UINT64_MAX + 1 - MOD : UINT64_MAX +
                                                                                                       1 +
                                                                                                       tmp_sum_bc +
                                                                                                       d_minus_ab;
                res = (res >= MOD) ? (res - MOD) : res;
            } else {
                res = tmp_sum_bc + d_minus_ab + UINT64_MAX + 1 - MOD;
            }
        } else {
            if (minus_flag) {
                res = (tmp_sum_bc >= (-d_minus_ab)) ? tmp_sum_bc + d_minus_ab : MOD + tmp_sum_bc + d_minus_ab;
            } else {
                res = ((tmp_sum_bc + d_minus_ab) >= MOD) ? tmp_sum_bc + d_minus_ab - MOD : tmp_sum_bc + d_minus_ab;
            }
        }
        return res;

    }

//https://1drv.ms/o/s!AvAq0B6vaN_Tj02p-JtAH0A9ICGw?e=ITjEAn
    Ntt64 fastmm_opt(Ntt64 x, Ntt64 y) {
        uint32_t x0 = (uint32_t) x;
        uint32_t x1 = (uint32_t) (x >> 32);
        uint32_t y0 = (uint32_t) y;
        uint32_t y1 = (uint32_t) (y >> 32);

        uint64_t x0y0 = (uint64_t) x0 * (uint64_t) y0;
        uint64_t x0y1 = (uint64_t) x0 * (uint64_t) y1;
        uint64_t x1y0 = (uint64_t) x1 * (uint64_t) y0;
        uint64_t x1y1 = (uint64_t) x1 * (uint64_t) y1;


        // z = x*y, z<127:0>
        // a = z<127:96>, b = z<95:64>, c = z<63:32>, d = z<31:0>
        auto d = (uint32_t) x0y0;
        uint64_t pp1 = (x0y0 >> 32) + (uint32_t) x1y0 + (uint32_t) x0y1;
        auto c = (uint32_t) pp1;
        uint64_t pp2 = (pp1 >> 32) + (x0y1 >> 32) + (x1y0 >> 32) + (uint32_t) x1y1;
        auto b = uint32_t(pp2);
        uint32_t a = (x1y1 >> 32) + (pp2 >> 32);
        uint64_t sum_bc = (uint64_t) b + (uint64_t) c;
        bool bc_of = (bool) (sum_bc >> 32);
        sum_bc = (uint64_t) (sum_bc << 32);
        uint64_t sum = sum_bc + (uint64_t) d;
        uint64_t minus = (uint64_t) a + (uint64_t) b;
        uint64_t res = 0;
        if (bc_of) {
            res = ((sum + UINT32_MAX) >= minus) ? sum - minus + UINT32_MAX : sum - minus + UINT32_MAX + MOD;
//        res =  sum - minus + UINT32_MAX;
//        res = (res>MOD)? res - MOD : res;
        } else {
            res = (sum >= minus) ? (((sum - minus) >= MOD) ? sum - minus - MOD : sum - minus) : MOD + sum - minus;
        }

        return res;

    }

    void genTW_ROM(TwRom &tw_rom) {
        auto w_n = (tw_rom.N) >> 1;
        auto phi_n = tw_rom.N;
        auto w_q = Ntt64((MOD - 1) / (w_n << 1));
        auto phi_q = Ntt64((MOD - 1) / (phi_n << 1));
        Ntt64 temp = 0;
        for (int i = 0; i < w_n; i++) {
            temp = POW(PRIM_ROOT, Ntt64(i * w_q), MOD);
            tw_rom.w_rom[i] = temp;
            tw_rom.inv_w_rom[i] = modINV(temp);
        }
        for (int j = 0; j < phi_n; j++) {
            temp = POW(PRIM_ROOT, Ntt64(j * phi_q), MOD);
            tw_rom.phi_rom[j] = temp;
            tw_rom.inv_phi_rom[j] = modINV(temp);
        }
    }

    void genNWCparam(TwParam &nwc_tw, const int n, const TwRom &tw_rom, const std::string &str) {
        int lvl = calLogBase2(n);
        auto w_n = n >> 1;
        auto phi_n = n;
//    Ntt64 w_q = Ntt64((MOD - 1)/(w_n<<1));
//    Ntt64 phi_q = Ntt64((MOD - 1)/(phi_n<<1));
        Ntt64 tw_temp = 0, phi_temp = 0, nwc_temp = 0;
        auto scale = 0, phi_size = 0, tw_size = 0, phi_probe = 0;
//    auto debug_tw = 0;
        if (str == "NWC-DIT-NR-NTT") {
            for (int i = 0; i < lvl; i++) {
                tw_size = n >> (lvl - i);
                phi_probe = n >> (i + 1);
                scale = w_n >> i;
                for (int j = 0; j < tw_size; j++) {
//                debug_tw = j*scale;
                    tw_temp = tw_rom.w_rom[j * scale];
                    phi_temp = tw_rom.phi_rom[phi_probe];
                    nwc_temp = modMULT(tw_temp, phi_temp);
                    nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);

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
//                debug_tw = j*scale;
                    tw_temp = tw_rom.inv_w_rom[j * scale];
                    phi_temp = tw_rom.inv_phi_rom[phi_probe];
                    nwc_temp = modMULT(tw_temp, phi_temp);
                    nwc_tw.tw_factor[i].push_back(nwc_temp);
//                nwc_tw.tw_factor[i].push_back(tw_temp);
                }
                bitRev(nwc_tw.tw_factor[i]);
            }
        }

    }

// Function to perform Number Theoretic Transform (NTT)
    void applyNtt(LagrangePolynomial &out, const IntPolynomial &in) {
        auto N = in.N;
        LagrangePolynomial format_input(N);
        for (int i = 0; i < N; i++) {
            if (in.coeffs[i] >= 0) {
                format_input.coeffs[i] = Ntt64(in.coeffs[i]);
            } else {
                format_input.coeffs[i] = Ntt64(in.coeffs[i] + MOD);
            }
        }
        DITNRLaPoly(out, format_input);
    }

    void applyIntt(IntPolynomial &out, const LagrangePolynomial &in) {
        auto N = in.N;
        LagrangePolynomial res(N);
        DIFRNLaPoly(res, in);
        int64_t temp_ntt;
        int64_t temp_poly;
        for (int i = 0; i < N; i++) {
            if (res.coeffs[i] >= HALF_MOD) {
                temp_ntt = int64_t(res.coeffs[i] - MOD);
            } else {
                temp_ntt = int64_t(res.coeffs[i]);
            }
            temp_poly = temp_ntt % TORUS_Q;
            if (temp_poly < TORUS_MIN) {
                out.coeffs[i] = int32_t(temp_poly + TORUS_Q);
            } else if (temp_poly > TORUS_MAX) {
                out.coeffs[i] = int32_t(temp_poly - TORUS_Q);
            } else {
                out.coeffs[i] = int32_t(temp_poly);
            }
        }
    }

// output_j = aj * bj mod p
    void modularMult(std::vector<uint64_t> &output, const std::vector<uint64_t> &coeffsA,
                     const std::vector<uint64_t> &coeffsB) {
        const auto N = output.size();
        for (auto j = 0; j < N; j++) {
            output[j] = fastmm_opt(coeffsA[j], coeffsB[j]);
        }
    }

// b += a * s mod p
    void
    modularAccumulate(std::vector<uint64_t> &res, const std::vector<uint64_t> &in1, const std::vector<uint64_t> &in2) {
        const auto N = res.size();
        Ntt64 tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (in1[j] == 0 || in2[j] == 0) {
                continue;
            }
            tmp = fastmm_opt(in1[j], in2[j]);
            res[j] = modAdd(res[j], tmp);
        }
    }

// res = aN * sN
    void calModularInnerProductNtt(LagrangePolynomial &res, const vector<LagrangePolynomial> &in1,
                                   const vector<LagrangePolynomial> &in2) {
        for (auto i = 0; i < in1.size(); i++) {
            modularAccumulate(res.coeffs, in1[i].coeffs, in2[i].coeffs);
        }
    }

    void
    calModularInnerProductNtt(LagrangePolynomial &res, const LagrangePolynomial &in1, const LagrangePolynomial &in2) {
        modularAccumulate(res.coeffs, in1.coeffs, in2.coeffs);
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