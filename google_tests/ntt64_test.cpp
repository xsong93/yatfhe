//
// Created by ic on 24-4-30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include <iostream>
#include <gmp.h>
#include "yautil/tool.h"
#include "yatfhe/ntt.h"
using namespace std;

TEST(ntt64_test, ntt64_test){
    int N = 4096;
    int depth = clog2(N);
    TW_PARAM nwc_tw(depth);
    TW_PARAM nwc_itw(depth);
    TW_ROM tw_rom(N);
    genTW_ROM(tw_rom);
    std::string str_ntt = "NWC-DIT-NR-NNT";
    std::string str_intt = "NWC-DIF-RN-INNT";
    genNWCparam(nwc_tw,N,tw_rom,str_ntt);
    genNWCparam(nwc_itw,N,tw_rom,str_intt);
    IntPolynomial a{N}, b{N}, res{N}, ref{N}, test{N};
    Ntt64Polynomial a_ntt{N}, b_ntt{N}, mul_ntt{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i - (N>>1);
        b.coeffs[i] = i + (1<<30);
    }
    polynomialMulNaive(ref,a,b);
    doNTT(a_ntt,a,nwc_tw);
    doNTT(b_ntt,b,nwc_tw);
    for (int i = 0; i < N; i++) {
        mul_ntt.coeffs[i] = modMULT(a_ntt.coeffs[i], b_ntt.coeffs[i]);
    }
    doINTT(res,mul_ntt,nwc_itw);
    for (int i = 0; i < N; i++) {
        EXPECT_EQ(res.coeffs[i], ref.coeffs[i]);
    }
    cout<<"breakpoint"<<endl;
    printBanner("NTT_NWC_64");
}

TEST(ntt64_test, single_test) {
    int N = 1024;
    int depth = clog2(N);
    TW_PARAM nwc_tw(depth);
    TW_PARAM nwc_itw(depth);
    TW_ROM tw_rom(N);
    genTW_ROM(tw_rom);
    std::string str_ntt = "NWC-DIT-NR-NNT";
    std::string str_intt = "NWC-DIF-RN-INNT";
    genNWCparam(nwc_tw,N,tw_rom,str_ntt);
    genNWCparam(nwc_itw,N,tw_rom,str_intt);
    Ntt64Polynomial a{N}, a_ntt{N}, res{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i;
    }
    DIT_NR(a_ntt, a, nwc_tw);
    DIF_RN(res, a_ntt, nwc_itw);
    cout<<"breakpoint"<<endl;
}

TEST(MODMULT_TEST, test1){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0;

    int N = 4096;
    int lvl = clog2(N);
    for (int i = 0 ; i < N; i++ ) {
        x += ((MOD64+1)>>lvl);
        y = x + 1;
        res1 = modMULT(x,y);
        res2 = modMul(x,y);
        EXPECT_EQ(res1, res2);
    }
    printBanner("modmult test");

}