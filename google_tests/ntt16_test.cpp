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
    auto add2 = modADD64(a, b);
    auto sub = modSUB16(a, b);
    auto sub2 = modSUB64(a, b);
    auto mul = modMULT16(a, b);
    auto mul2 = modMULT64(a, b);
    auto inv = modINV16(a);
    auto inv2 = modINV64(a);
    printf("a+b: %d, a-b: %d, a*b: %d, 1/a: %d\n", add, sub, mul, inv);
    printf("a+b: %d, a-b: %d, a*b: %d, 1/a: %d\n", add2, sub2, mul2, inv2);
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