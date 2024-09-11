//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt_old.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

TEST(NttOldTest, NttOldBasicArithTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 1024;
    LagrangePolynomial a{N};
    LagrangePolynomial b{N};
    LagrangePolynomial tmpMul{N};
    LagrangePolynomial tmpAdd{N};
    LagrangePolynomial tmpSub{N};

    TorusPolynomial poly0{N};
    TorusPolynomial poly2{N};
    TorusPolynomial resMul{N};
    TorusPolynomial resAdd{N};
    TorusPolynomial resSub{N};
    TorusPolynomial navMul{N};
    TorusPolynomial navAdd{N};
    TorusPolynomial navSub{N};
    int t = 10;
    while (t-- > 0) {
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = genIntUniformDist(1 << 25, 1<< 30);
            poly2.coeffs[j] = genIntUniformDist(INT32_MIN, INT32_MAX);
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly2.coeffs, "poly2");

        COUNT_TIME("NTT_MULT", {
            applyNttOld(a, poly0);
            applyNttOld(b, poly2);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = modMulOld(a.coeffs[i], b.coeffs[i]);
            }
            applyInttOld(resMul, tmpMul);
        })
        COUNT_TIME("NAIVE_MULT",
            polynomialMulNaive(navMul, poly0, poly2);)
            for (int i = 0; i < a.N; i++) {
                tmpAdd.coeffs[i] = modAddOld(a.coeffs[i], b.coeffs[i]);
                tmpSub.coeffs[i] = modSubOld(a.coeffs[i], b.coeffs[i]);
        }
        applyInttOld(resAdd, tmpAdd);
        applyInttOld(resSub, tmpSub);

        polynomialAdd(navAdd, poly0, poly2);
        polynomialSub(navSub, poly0, poly2);

        for (int i = 0; i < navMul.N; i++) {
            EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
//            EXPECT_EQ(resAdd.coeffs[i], navAdd.coeffs[i]);
//            EXPECT_EQ(resSub.coeffs[i], navSub.coeffs[i]);
        }
    }
    printBanner("NttBasicArithTest");
}
