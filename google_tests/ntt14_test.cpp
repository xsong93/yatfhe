//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include <gmp.h>
#include "yatfhe/ntt.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

using namespace NttNative14;

TEST(NTT14, MOD_MULT) {
    const int N = 512;
    initGlobalParamsNtt(N);
    Ntt14 a = 12280;
    Ntt14 b = 11111;
    auto add = modADD(a, b);
    auto sub = modSUB(a, b);
    auto mul = modMULT(a, b);
    auto inv = modINV(a);
    auto pow = POW(a, 3);
    printf("a+b: %d, a-b: %d, a*b: %d, inv: %d, pow: %d\n", add, sub, mul, inv, pow);
    ASSERT_EQ(add, 11102);
    ASSERT_EQ(sub, 1169);
    ASSERT_EQ(mul, 10602);
    ASSERT_EQ(inv, 9558);
    ASSERT_EQ(pow, ((uint32_t (a * a) % MOD) * a) % MOD);

    printBanner("NTT14.MOD_MULT");
}

TEST(NTT14, NTT_INTT) {
    const int N = 1024;
    initGlobalParamsNtt(N);
    Ntt14Polynomial resNtt{N};
    IntPolynomial a1{N};
    IntPolynomial resIntt{N};
    for (auto i = 0; i < N; i++) {
        a1.coeffs[i] = genIntUniformDist(CHAR_MIN, CHAR_MAX);
//        a1.coeffs[i] = 65536;
    }

    applyNtt(resNtt, a1);
    applyIntt(resIntt, resNtt);
    printArray(resNtt.coeffs, "resNtt");
    printArray(a1.coeffs, "orig");
    printArray(resIntt.coeffs, "intt");
    for (auto i = 0; i < a1.N; i++) {
        if (a1.coeffs[i] != resIntt.coeffs[i]) {
            printf("NE at index %d\n",i);
        }
        ASSERT_EQ(a1.coeffs[i], resIntt.coeffs[i]);
    }
    printBanner("NTT14.NTT_INTT");
}