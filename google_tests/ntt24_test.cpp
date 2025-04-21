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

using namespace NttNative24;

TEST(Ntt24Test, PowInvTest) {
    std::vector<Ntt24> a(0);
    std::vector<Ntt24> b(0);
    std::vector<Ntt24> c(0);
    for (int i = 1; i < 100; i++) {
        a.push_back(POW(i, i));
        b.push_back(modINV(a[i - 1]));
        c.push_back(modInverse(a[i-1], MOD));
    }
    printArray(a, "a");
    printArray(b, "b");
    printArray(c, "c");
    printBanner("PowInvTest");
}

TEST(Ntt24Test, NttIntt24Test) {
    const int N = 512;
    initGlobalParamsNtt(N);
    Ntt24Polynomial resNtt{N};
    Int8Polynomial a1{N};
    Int8Polynomial resIntt{N};
    int q = 251;
    for (auto i = 0; i < N; i++) {
        a1.coeffs[i] = static_cast<int8_t>(genIntUniformDist(-q/2, (q - 1)/2));
//        a1.coeffs[i] = 65536;
    }

    applyNtt(resNtt, a1);
    applyIntt(resIntt, resNtt, q);
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
    initGlobalParamsNtt(N);
    Ntt24Polynomial a{N};
    Ntt24Polynomial b{N};
    Ntt24Polynomial tmpMul{N};

    Int8Polynomial poly0{N};
    Int8Polynomial poly2{N};
    Int8Polynomial resMul{N};
    Int8Polynomial navMul{N};
    int q = 251;
    int t = 0;
    while (t++ < 10) {
        cout << "cycle: " << t << endl;
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = static_cast<int8_t>(genIntUniformDist(-q/2, (q - 1)/2));
            poly2.coeffs[j] = static_cast<int8_t>(genIntUniformDist(-q/2, (q - 1)/2));
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly2.coeffs, "poly2");

        COUNT_TIME("NTT_MULT", {
            applyNtt(a, poly0);
            applyNtt(b, poly2);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = modMULT(a.coeffs[i], b.coeffs[i]);
            }
            applyIntt(resMul, tmpMul, q);
        })
        COUNT_TIME("NAIVE_MULT",
                   multInt8Polynomial(navMul, poly0, poly2, q);)

        printArray(resMul.coeffs, "resMul");
        printArray(navMul.coeffs, "navMul");

        for (int i = 0; i < navMul.N; i++) {
            EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        }
    }
    printBanner("Ntt24BasicArithTest");
}