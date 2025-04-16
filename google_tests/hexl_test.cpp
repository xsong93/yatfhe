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
#include "yatfhe/ntt.h"
#include <gtest/gtest.h>

TEST(HEXL_TEST, NTT_INTT) {
    YatfheParameters param {};
    param.N = 1024;
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    auto N = param.N;

    printHexlParams();
    TorusPolynomial in{N};
    LagrangePolynomial nttHexl{N};
    LagrangePolynomial ntt{N};
    TorusPolynomial outHexl{N};
    TorusPolynomial out{N};
    for (auto i = 0; i < N; i++) {
        in.coeffs[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    }

    COUNT_TIME("HEXL", applyNttHexl(nttHexl, in);)
    COUNT_TIME("HEXL", applyInttHexl(outHexl, nttHexl);)
    COUNT_TIME("n32", applyNtt(ntt, in);)
    COUNT_TIME("n32", applyIntt(out, ntt);)

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
//            poly2.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
//            poly3.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
//            poly4.coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly1.coeffs, "poly1");
//        printArray(poly2.coeffs, "poly2");
//        printArray(poly3.coeffs, "poly3");
//        printArray(poly4.coeffs, "poly4");

        COUNT_TIME("HEXL_MULT", {
            applyNttHexl(a, poly0);
            applyNttHexl(b, poly1);
//            applyNtt32(c, poly2);
//            applyNtt32(d, poly3);
//            applyNtt32(e, poly4);
            EltwiseMultMod(tmpMul1.coeffs.data(), a.coeffs.data(), b.coeffs.data(), N, p.qNtt, 1);
            applyInttHexl(resMul1, tmpMul1);
        })
        COUNT_TIME("NTT_MULT", {
            applyNtt(a, poly0);
            applyNtt(b, poly1);
//            applyNtt32(c, poly2);
//            applyNtt32(d, poly3);
//            applyNtt32(e, poly4);
            calModularInnerProductNtt(tmpMul2, a, b);
            applyIntt(resMul2, tmpMul2);
        })
        COUNT_TIME("NAIVE_MULT", {
            TorusPolynomial tmp{N};
            TorusPolynomial tmp1{N};
            polynomialMulNaiveModQ(navMul, poly0, poly1, TORUS_Q);
//            polynomialMulNaiveModQ(tmp1, tmp, poly2, TORUS_Q);
//            polynomialMulNaiveModQ(tmp, tmp1, poly3, TORUS_Q);
//            polynomialMulNaiveModQ(navMul, tmp, poly4, TORUS_Q);
        })

        printArray(resMul1.coeffs, "resMul1");
        printArray(resMul2.coeffs, "resMul2");
        printArray(navMul.coeffs, "navMul");

        EXPECT_EQ(resMul1.coeffs, navMul.coeffs);
    }
    printBanner("Ntt32BasicArithTest");
}