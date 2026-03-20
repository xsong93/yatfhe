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

TEST(Ntt14Test, ModMultTest) {
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
}

TEST(Ntt14Test, NttIntt14Test) {
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
    printBanner("NttIntt14Test");
}

TEST(Ntt14Test, Ntt14BasicArithTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 512;
    initGlobalParamsNtt(N);
    Ntt14Polynomial a{N};
    Ntt14Polynomial b{N};
    Ntt14Polynomial tmpMul{N};

    Int8Polynomial poly0{N};
    Int8Polynomial poly2{N};
    Int8Polynomial resMul{N};
    Int8Polynomial navMul{N};
    int t = 1;
    while (t-- > 0) {
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = genIntUniformDist((1 << 4), (1 << 4));
            poly2.coeffs[j] = genIntUniformDist((1 << 4), (1 << 4));
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly2.coeffs, "poly2");

        COUNT_TIME("NTT_MULT", {
            applyNttPoly8(a, poly0);
            applyNttPoly8(b, poly2);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = modMULT(a.coeffs[i], b.coeffs[i]);
            }
            applyInttPoly8(resMul, tmpMul);
        })
        COUNT_TIME("NAIVE_MULT",
                   multInt8Polynomial(navMul, poly0, poly2, 1 << 8);)

        printArray(resMul.coeffs, "resMul");
        printArray(navMul.coeffs, "navMul");

        for (int i = 0; i < navMul.N; i++) {
            EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        }
    }
    printBanner("NttSamePoly");
}