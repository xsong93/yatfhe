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
    Ntt32_TW TW(512);
    Ntt32_iTW iTW(512);
    Ntt32 test_res = POW(7,10000,MOD);
    std::cout<<test_res<<std::endl;
    std::cout<<clog2(1024)<<std::endl;
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
    NttPolynomial a{1024},b {1024};
    for (int i = 0; i < 1024; i++) {
        a.coeffs[i] = i;
    }
    std::cout<<"a:N = "<<a.N<<"; b:N = "<<b.N<<std::endl;
    doNTT32(b,a,TW);
    printNttPoly(b);
}