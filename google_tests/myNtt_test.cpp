//
// Created by ic on 24-4-11.
//
#include "gtest/gtest.h"
#include "yatfhe/myNtt.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include <gmp.h>
#include "yautil/tool.h"

TEST(MYNTT_TEST,test_init) {
    int N = 1024;
    Ntt32_TW TW(N>>1);
    Ntt32_iTW iTW(N>>1);
    std::cout<<MOD<<std::endl;
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
    NttPolynomial a{N},b {N}, res{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i;
    }
    std::cout<<"a:N = "<<a.N<<"; b:N = "<<b.N<<std::endl;
    doNTT32(b,a,TW);
    printNttPoly(b);
    doINTT32(res, b, iTW);
    printNttPoly(res);
}

TEST(MYNTT_TEST, modSUB_test){
    Ntt32 a = 0, b = 1;
    std::cout<<modSUB(a,b)<<std::endl;
}
TEST(MYNTT_TEST, modADD_test){
    Ntt32 a = 4293918719, b = 74203740;
    std::cout<<modADD(a,b)<<std::endl;
}

TEST(MYNTT_TEST, mod_inv_test) {
    int N = 4;
    Ntt32_TW TW(N>>1);
    Ntt32_iTW iTW(N>>1);
    std::cout<<MOD<<std::endl;
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
}

TEST(MYNTT_TEST, modMULT_test) {
    Ntt32 a = 4293918719, b = 37101870;
    std::cout<< modMULT(a,b)<<std::endl;
}