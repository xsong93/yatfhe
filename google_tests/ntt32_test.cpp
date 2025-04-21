//
// Created by Xintong Song on 2032/1/30.
//
#include "gtest/gtest.h"
#include <gmp.h>
#include "yatfhe/ntt32.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

using namespace NttNative32;

TEST(Ntt32Test, PowInvTest) {
    std::vector<Ntt32> a(0);
    std::vector<Ntt32> b(0);
    std::vector<Ntt32> c(0);
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

TEST(Ntt32Test, NttIntt32Test) {
    YatfheParameters p{};
    p.q = Q_32P;
    initYatfhe(p);

    auto N = p.N;
    Ntt32Polynomial resNtt{N};
    TorusPolynomial a1{N};
    TorusPolynomial resIntt{N};
    for (auto i = 0; i < N; i++) {
        a1.coeffs[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
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
    printBanner("NttIntt32Test");
}

TEST(Ntt32Test, Ntt32BasicArithTest) {
    YatfheParameters p{};
    p.q = Q_32P;
    p.N = 1024;
    initYatfhe(p);

    auto N = p.N;
    Ntt32Polynomial a{N};
    Ntt32Polynomial b{N};
    Ntt32Polynomial c{N};
    Ntt32Polynomial d{N};
    Ntt32Polynomial e{N};
    Ntt32Polynomial tmpMul{N};

    TorusPolynomial poly0{N};
    TorusPolynomial poly1{N};
    TorusPolynomial poly2{N};
    TorusPolynomial poly3{N};
    TorusPolynomial poly4{N};
    TorusPolynomial resMul{N};
    TorusPolynomial navMul{N};
    int t = 0;
    while (t++ < 1) {
        cout << "cycle: " << t << endl;
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
            poly1.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
            poly2.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
            poly3.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
            poly4.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly1.coeffs, "poly1");
        printArray(poly2.coeffs, "poly2");
        printArray(poly3.coeffs, "poly3");
        printArray(poly4.coeffs, "poly4");

        COUNT_TIME("NTT_MULT", {
            applyNtt(a, poly0);
            applyNtt(b, poly1);
            applyNtt(c, poly2);
            applyNtt(d, poly3);
            applyNtt(e, poly4);
            for (int i = 0; i < a.N; i++) {
                auto tmp = modMULT(a.coeffs[i], b.coeffs[i]);
                tmp = modMULT(tmp, c.coeffs[i]);
                tmp = modMULT(tmp, d.coeffs[i]);
                tmpMul.coeffs[i] = modMULT(tmp, e.coeffs[i]);
            }
            applyIntt(resMul, tmpMul);
        })
        COUNT_TIME("NAIVE_MULT", {
            TorusPolynomial tmp{N};
            TorusPolynomial tmp1{N};
            multIntPolynomialModQ(tmp, poly0, poly1, TORUS_Q);
            multIntPolynomialModQ(tmp1, tmp, poly2, TORUS_Q);
            multIntPolynomialModQ(tmp, tmp1, poly3, TORUS_Q);
            multIntPolynomialModQ(navMul, tmp, poly4, TORUS_Q);
        })

        printArray(resMul.coeffs, "resMul");
        printArray(navMul.coeffs, "navMul");

        EXPECT_EQ(resMul.coeffs, navMul.coeffs);
    }
    printBanner("Ntt32BasicArithTest");
}