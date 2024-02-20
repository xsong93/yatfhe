//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yautil/time_counter.h"

TEST(NttTest, X) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 1024;
    LagrangePolynomial a(N);
    LagrangePolynomial b(N);
    LagrangePolynomial tmpMul(N);
    LagrangePolynomial tmpAdd(N);
    LagrangePolynomial tmpSub(N);

    IntPolynomial poly1(N);
    TorusPolynomial poly2(N);
    TorusPolynomial resMul(N);
    TorusPolynomial resAdd(N);
    TorusPolynomial resSub(N);
    TorusPolynomial navMul(N);
    TorusPolynomial navAdd(N);
    TorusPolynomial navSub(N);

    for (int i = 0; i < a.N; i++) {
        poly1.coeffs[i] = i;
        poly2.coeffs[i] = i;
    }
    COUNT_TIME("NTT_MULT",
                    applyNtt(a, poly1);
                    applyNtt(b, poly2);
                    for (int i = 0; i < a.N; i++) {
                       tmpMul.coeffs[i] = modMul(a.coeffs[i], b.coeffs[i]);
                    }
                    applyIntt(resMul, tmpMul);)
    COUNT_TIME("NAIVE_MULT",
                    polynomialMulNaive(navMul, poly1, poly2);)

    for (int i = 0; i < a.N; i++) {
        tmpAdd.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
        tmpSub.coeffs[i] = modSub(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(resAdd, tmpAdd);
    applyIntt(resSub, tmpSub);

    polynomialAdd(navAdd, poly1, poly2);
    polynomialSub(navSub, poly1, poly2);
    
    for (int i = 0; i < navMul.N; i++) {
        EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        EXPECT_EQ(resAdd.coeffs[i], navAdd.coeffs[i]);
        EXPECT_EQ(resSub.coeffs[i], navSub.coeffs[i]);
    }
    std::cout << ">>>>>>>> NTT test passed! <<<<<<<<" << std::endl;
}