//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/trgsw.h"
#include "yautil/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

TEST(NttAddConstantTest, NttAddConstantTest) {
    const int N = 1024;
    LagrangePolynomial a(N);
    LagrangePolynomial b(N);
    LagrangePolynomial resNtt(N);
    IntPolynomial poly(N);
    IntPolynomial c(N);
    IntPolynomial res(N);
    for (int i = 0; i < a.N; i++) {
        poly.coeffs[i] = i;
    }
    c.coeffs[1] = 77;
    applyNtt(a, poly);
    applyNtt(b, c);
    printArray(c.coeffs, "cOri");
    printArray(b.coeffs, "bNtt");
    for (int i = 0; i < a.N; i++) {
        resNtt.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(res, resNtt);
    printArray(res.coeffs, "res");
}

TEST(NttSamePolyTest, NttSamePolyTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 1024;
    LagrangePolynomial a(N);
    LagrangePolynomial b(N);
    LagrangePolynomial tmpMul(N);
    LagrangePolynomial tmpAdd(N);
    LagrangePolynomial tmpSub(N);

    TorusPolynomial poly0(N);
    TorusPolynomial poly2(N);
    TorusPolynomial resMul(N);
    TorusPolynomial resAdd(N);
    TorusPolynomial resSub(N);
    TorusPolynomial navMul(N);
    TorusPolynomial navAdd(N);
    TorusPolynomial navSub(N);

    initCoeffsViaUniformDistribution(poly0.coeffs, N);
    initCoeffsViaUniformDistribution(poly2.coeffs, N);
    printArray(poly0.coeffs, "poly0");
    printArray(poly2.coeffs, "poly2");

    COUNT_TIME("NTT_MULT",
               applyNtt(a, poly0);
                       applyNtt(b, poly2);
                       for (int i = 0; i < a.N; i++) {
                           tmpMul.coeffs[i] = modMul(a.coeffs[i], b.coeffs[i]);
                       }
                       applyIntt(resMul, tmpMul);)
    COUNT_TIME("NAIVE_MULT",
               polynomialMulNaive(navMul, poly0, poly2);)

    for (int i = 0; i < a.N; i++) {
        tmpAdd.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
        tmpSub.coeffs[i] = modSub(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(resAdd, tmpAdd);
    applyIntt(resSub, tmpSub);

    polynomialAdd(navAdd, poly0, poly2);
    polynomialSub(navSub, poly0, poly2);

    for (int i = 0; i < navMul.N; i++) {
        EXPECT_NEAR(resMul.coeffs[i], navMul.coeffs[i], doubleToTorus32(0.01));
        EXPECT_EQ(resAdd.coeffs[i], navAdd.coeffs[i]);
        EXPECT_EQ(resSub.coeffs[i], navSub.coeffs[i]);
    }
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> NttSamePolyTest test passed! <<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}

TEST(NttDiffPolyTest, NttDiffPolyTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 1024;
    const int mSize = 2 * N;

    LagrangePolynomial a(N);
    LagrangePolynomial b(N);
    LagrangePolynomial tmpMul(N);

    TorusPolynomial poly0(N);
    IntPolynomial poly01(N);
    TorusPolynomial resMul(N);
    TorusPolynomial navMul(N);

    uniform_int_distribution<int> distribution(0, 1);
    for (int j = 0; j < N; j++) {
        poly01.coeffs[j] = distribution(rng);
    }
    printArray(poly01.coeffs, "poly01");

    initCoeffsViaUniformDistribution(poly0.coeffs, N);
    printArray(poly0.coeffs, "poly0");

    COUNT_TIME("NTT_MULT",
               applyNtt(a, poly0);
               applyNtt(b, poly01);
               for (int i = 0; i < a.N; i++) {
                   tmpMul.coeffs[i] = modMul(a.coeffs[i], b.coeffs[i]);
               }
               applyIntt(resMul, tmpMul);)
    COUNT_TIME("NAIVE_MULT",
               polynomialMulNaive(navMul, poly0, poly01);)

    for (int i = 0; i < navMul.N; i++) {
        EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
    }
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> NTT test passed! <<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}