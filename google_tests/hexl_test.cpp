//
// Created by Xintong  on 25-4-16.
//
#include <iostream>
#include "yatfhe/ntt_hexl.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/initializer.h"
#include "yatfhe/ntt64.h"
#include <gtest/gtest.h>

TEST(HEXL_TEST, NTT_INTT) {
    YatfheParameters param {};
    param.N = 1024;
    initYatfhe(param);

    auto N = param.N;
    TorusPolynomial in{N};
    NttPolynomial nttHexl{N};
    NttPolynomial ntt{N};
    TorusPolynomial outHexl{N};
    TorusPolynomial out{N};
    for (auto i = 0; i < N; i++) {
        in.coeffs[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    }

    COUNT_TIME("HEXL_NTT", NttHexl::applyNtt(nttHexl, in);)
    COUNT_TIME("HEXL_INTT", NttHexl::applyIntt(outHexl, nttHexl);)
    COUNT_TIME("NATIVE_NTT", NttNative64::applyNtt(ntt, in);)
    COUNT_TIME("NATIVE_INTT", NttNative64::applyIntt(out, ntt);)

    printArray(in.coeffs, "in");
    printArray(outHexl.coeffs, "outHexl");
    printArray(out.coeffs, "out");
    ASSERT_EQ(out.coeffs, outHexl.coeffs);
}

TEST(HEXL_TEST, POLY_MULT) {
    YatfheParameters p{};
    p.N = 1024;
    initYatfhe(p);

    auto N = p.N;
    NttPolynomial a{N};
    NttPolynomial b{N};
    NttPolynomial c{N};
    NttPolynomial d{N};
    NttPolynomial e{N};
    NttPolynomial tmpMul1{N};
    NttPolynomial tmpMul2{N};

    TorusPolynomial poly0{N};
    TorusPolynomial poly1{N};
    TorusPolynomial poly2{N};
    TorusPolynomial poly3{N};
    TorusPolynomial poly4{N};
    TorusPolynomial resMul1{N};
    TorusPolynomial resMul2{N};
    TorusPolynomial navMul{N};
    int t = 0;
    while (t++ < 1) {
        cout << "cycle: " << t << endl;
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
            poly1.coeffs[j] = genIntUniformDist(0, 1);
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly1.coeffs, "poly1");

        COUNT_TIME("HEXL_MULT", {
            NttHexl::applyNtt(a, poly0);
            NttHexl::applyNtt(b, poly1);
            EltwiseMultMod(tmpMul1.coeffs.data(), a.coeffs.data(), b.coeffs.data(), N, p.qNtt, 1);
            NttHexl::applyIntt(resMul1, tmpMul1);
        })
        COUNT_TIME("NTT_MULT", {
            NttNative64::applyNtt(a, poly0);
            NttNative64::applyNtt(b, poly1);
            NttNative64::calModularInnerProductNtt(tmpMul2, a, b);
            NttNative64::applyIntt(resMul2, tmpMul2);
        })
        COUNT_TIME("NAIVE_MULT", {
            TorusPolynomial tmp{N};
            TorusPolynomial tmp1{N};
            multIntPolynomialModQ(navMul, poly0, poly1, TORUS_Q);
        })

        printArray(resMul1.coeffs, "resMul1");
        printArray(resMul2.coeffs, "resMul2");
        printArray(navMul.coeffs, "navMul");

        EXPECT_EQ(resMul1.coeffs, navMul.coeffs);
    }
    printBanner("Ntt32BasicArithTest");
}

TEST(HEXL_TEST, NTT_ROT) {
    YatfheParameters param {};
    // param.N = 4;
    // param.qNtt = 7681;
    initYatfhe(param);

    auto N = param.N;

    TorusPolynomial in{N};
    TorusPolynomial in2{N};
    TorusPolynomial ref{N};
    NttPolynomial nttHexl{N}, nttHexl2{N};
    NttPolynomial tmp{N};
    TorusPolynomial res{N};
    for (auto i = 0; i < N; i++) {
        in.coeffs[i] = genIntUniformDist(-4, 4);
        in2.coeffs[i] = -in.coeffs[i];
    }
    int r = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    cout << "r:" << r << endl;
    COUNT_TIME("torusPolynomialRotate", rotateTorusPolynomial(ref, r, in));
    COUNT_TIME("HEXL", NttHexl::applyNtt(nttHexl, in);)
    COUNT_TIME("HEXL", NttHexl::applyNtt(nttHexl2, in2);)
    printArray(nttHexl.coeffs, "in");
    printArray(nttHexl2.coeffs, "in2");
    COUNT_TIME("NTT_ROT", rotateNttPolynomial(tmp, nttHexl, r);)
    COUNT_TIME("HEXL", NttHexl::applyIntt(res, tmp);)
    printArray(in.coeffs, "in");
    printArray(ref.coeffs, "ref");
    printArray(res.coeffs, "res");
    ASSERT_EQ(ref.coeffs, res.coeffs);
}
