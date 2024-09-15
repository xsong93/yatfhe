//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include <gmp.h>
#include "yatfhe/ntt24.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

TEST(Ntt24Test, PowInvTest) {
    std::vector<Ntt24> a(0);
    std::vector<Ntt24> b(0);
    std::vector<Ntt24> c(0);
    for (int i = 1; i < 100; i++) {
        a.push_back(POW24(i, i));
        b.push_back(modINV24(a[i-1]));
        c.push_back(modInverse(a[i-1], MOD24));
    }
    printArray(a, "a");
    printArray(b, "b");
    printArray(c, "c");
    printBanner("PowInvTest");
}

TEST(Ntt24Test, NttIntt24Test) {
    const int N = 512;
    initGlobalParamsNtt24(N);
    Ntt24Polynomial resNtt{N};
    Int8Polynomial a1{N};
    Int8Polynomial resIntt{N};
    for (auto i = 0; i < N; i++) {
        a1.coeffs[i] = static_cast<int8_t>(genIntUniformDist(CHAR_MIN, CHAR_MAX));
//        a1.coeffs[i] = 65536;
    }

    applyNtt24(resNtt, a1);
    applyIntt24(resIntt, resNtt);
    printArray(resNtt.coeffs, "resNtt");
    printArray(a1.coeffs, "orig");
    printArray(resIntt.coeffs, "intt");
    for (auto i = 0; i < a1.N; i++) {
        if (a1.coeffs[i] != resIntt.coeffs[i]) {
            printf("NE at index %d\n",i);
        }
        ASSERT_EQ(a1.coeffs[i], resIntt.coeffs[i]);
    }
    printBanner("NttIntt24Test");
}

TEST(Ntt24Test, Ntt24BasicArithTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 512;
    initGlobalParamsNtt24(N);
    Ntt24Polynomial a{N};
    Ntt24Polynomial b{N};
    Ntt24Polynomial tmpMul{N};

    Int8Polynomial poly0{N};
    Int8Polynomial poly2{N};
    Int8Polynomial resMul{N};
    Int8Polynomial navMul{N};
    int t = 0;
    while (t++ < 10) {
        cout << "cycle: " << t << endl;
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = static_cast<int8_t>(genIntUniformDist(CHAR_MIN, CHAR_MAX));
            poly2.coeffs[j] = static_cast<int8_t>(genIntUniformDist(CHAR_MIN, CHAR_MAX));
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly2.coeffs, "poly2");

        COUNT_TIME("NTT_MULT", {
            applyNtt24(a, poly0);
            applyNtt24(b, poly2);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = modMULT24(a.coeffs[i], b.coeffs[i]);
            }
            applyIntt24(resMul, tmpMul);
        })
        COUNT_TIME("NAIVE_MULT",
                   polynomialMulNaiveModQ8(navMul, poly0, poly2, 1 << 8);)

        printArray(resMul.coeffs, "resMul");
        printArray(navMul.coeffs, "navMul");

        for (int i = 0; i < navMul.N; i++) {
            EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        }
    }
    printBanner("Ntt24BasicArithTest");
}