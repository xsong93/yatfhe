//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

TEST(Ntt16Test, ModMultTest) {
    const int N = 512;
    initGlobalParamsNtt16(N);
    Ntt16 a = 59923;
    Ntt16 b = 65535;
    auto add = modADD16(a, b);
    auto sub = modSUB16(a, b);
    auto mul = modMULT16(a, b);
    auto inv = modINV16(a);
    auto pow = POW16(a, 3);
    printf("a+b: %d, a-b: %d, a*b: %d, inv: %d, pow: %d\n", add, sub, mul, inv, pow);
    ASSERT_EQ(add, 59921);
    ASSERT_EQ(sub, 59925);
    ASSERT_EQ(mul, 11228);
    ASSERT_EQ(inv, 36831);
    ASSERT_EQ(pow, ((uint32_t (a * a) % MOD16) * a) %MOD16);
}

TEST(Ntt16Test, NttIntt16Test) {
    const int N = 512;
    initGlobalParamsNtt16(N);
    Ntt16Polynomial resNtt{N};
    IntPolynomial a1{N};
    IntPolynomial resIntt{N};
    for (auto i = 0; i < N; i++) {
        a1.coeffs[i] = genIntUniformDist(CHAR_MIN, CHAR_MAX);
    }

    applyNtt16(resNtt, a1);
    applyIntt16(resIntt, resNtt);
    printArray(resNtt.coeffs, "resNtt");
    printArray(a1.coeffs, "orig");
    printArray(resIntt.coeffs, "intt");
    for (auto i = 0; i < a1.N; i++) {
        if (a1.coeffs[i] != resIntt.coeffs[i]) {
            printf("NE at index %d\n",i);
        }
        ASSERT_EQ(a1.coeffs[i], resIntt.coeffs[i]);
    }
    printBanner("NttIntt64Test");
}