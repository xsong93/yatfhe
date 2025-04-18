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
    yatfheInit(param);

    auto N = param.N;
    TorusPolynomial in{N};
    LagrangePolynomial nttHexl{N};
    LagrangePolynomial ntt{N};
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
    yatfheInit(p);

    auto N = p.N;
    LagrangePolynomial a{N};
    LagrangePolynomial b{N};
    LagrangePolynomial c{N};
    LagrangePolynomial d{N};
    LagrangePolynomial e{N};
    LagrangePolynomial tmpMul1{N};
    LagrangePolynomial tmpMul2{N};

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
            polynomialMulNaiveModQ(navMul, poly0, poly1, TORUS_Q);
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
    param.N = 1024;
//    param.qNtt = 7681;
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    auto N = param.N;
    auto q = param.qNtt;

    TorusPolynomial in{N};
    TorusPolynomial in1{N};
    TorusPolynomial in2{N};
    LagrangePolynomial nttHexl{N};
    LagrangePolynomial nttHexl1{N};
    LagrangePolynomial nttHexl2{N};
    LagrangePolynomial tmp{N};
    TorusPolynomial outHexl{N};
    for (auto i = 0; i < N; i++) {
        in.coeffs[i] = genIntUniformDist(-4, 4);
    }
    int r = 3;
    in1.coeffs[r] = 1;
    torusPolynomialRotate(in2, r, in);

    COUNT_TIME("HEXL", NttHexl::applyNtt(nttHexl, in);)
    COUNT_TIME("HEXL", NttHexl::applyNtt(nttHexl1, in1);)
    COUNT_TIME("HEXL", NttHexl::applyNtt(nttHexl2, in2);)

    printArray(nttHexl.coeffs, "ntt");
    printArray(nttHexl1.coeffs, "ntt1");
    printArray(nttHexl2.coeffs, "ntt2");

    COUNT_TIME("NTT_ROT",
        for (size_t i = 0; i < N; i++) {
            EltwiseMultMod(tmp.coeffs.data(), nttHexl.coeffs.data(), nttHexl1.coeffs.data(), N, q, 1);
        })

    COUNT_TIME("HEXL", NttHexl::applyIntt(outHexl, tmp);)

    printArray(in.coeffs, "in");
    printArray(outHexl.coeffs, "outHexl");
    ASSERT_EQ(in2.coeffs, outHexl.coeffs);
}
