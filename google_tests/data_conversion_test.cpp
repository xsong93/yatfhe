//
// Created by Xintong Song on 2024/2/29.
//
#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"


TEST(doubleTorusTest, doubleTorusTest) {
    int N = 100;
    DoublePolynomial doublePoly(N);
    TorusPolynomial tPoly(N);
    TorusPolynomial t2Poly(N);
    DoublePolynomial resPoly(N);

    double i = -0.5;
    for (auto& item : doublePoly.coeffs) {
        item = i;
        i += 0.01;
    }
    printArray(doublePoly.coeffs, "doublePoly");
    for (auto j = 0; j < N; j++) {
        tPoly.coeffs[j] = doubleToTorus32(doublePoly.coeffs[j]);
    }
    printArray(tPoly.coeffs, "torusPoly");
    for (auto j = 0; j < N; j++) {
        resPoly.coeffs[j] = torus32ToDouble(tPoly.coeffs[j]);
    }
    printArray(resPoly.coeffs, "resPoly");
    for (auto j = 0; j < N; j++) {
        EXPECT_NEAR(resPoly.coeffs[j], doublePoly.coeffs[j], 1e-3);
    }


    for (auto k = 0; k < 100; k++) {
        initCoeffsViaUniformDistribution(tPoly.coeffs);
//        printArray(tPoly.coeffs, "torusPoly");
        torusPolyToDoublePoly(resPoly, tPoly);
//        printArray(resPoly.coeffs, "resPoly");
        doublePolyToTorusPoly(t2Poly, resPoly);
//        printArray(t2Poly.coeffs, "t2Poly");
        for (auto j = 0; j < N; j++) {
            EXPECT_EQ(tPoly.coeffs[j], t2Poly.coeffs[j]);
        }
    }
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> Data(double/ Torus) conversion test passed! <<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}

TEST(intTorusTest, intTorusTest) {
    int N = 100;
    IntPolynomial intPoly(N);
    TorusPolynomial tPoly(N);
    IntPolynomial resPoly(N);
    TorusPolynomial resTPoly(N);
    int i = -50;
    for (auto& item : intPoly.coeffs) {
        item = i++;
    }
    printArray(intPoly.coeffs, "intPoly");
    intPolyToTorusPoly(tPoly, intPoly, N);
    printArray(tPoly.coeffs, "torusPoly");
    torusPolyToIntPoly(resPoly, tPoly, N);
    printArray(resPoly.coeffs, "resPoly");
    for (auto j = 0; j < N; j++) {
        EXPECT_EQ(resPoly.coeffs[j], intPoly.coeffs[j]);
    }

    for (auto k = 0; k < 100; k++) {
        initCoeffsViaUniformDistribution(tPoly.coeffs);
//        printArray(tPoly.coeffs, "torusPoly");
        torusPolyToIntPoly(resPoly, tPoly, 2 * N);
//        printArray(resPoly.coeffs, "resPoly");
        intPolyToTorusPoly(resTPoly, resPoly, 2 * N);
//        printArray(resTPoly.coeffs, "torusPoly");
        for (auto j = 0; j < N; j++) {
            EXPECT_NEAR(resTPoly.coeffs[j], tPoly.coeffs[j], doubleToTorus32(0.01));
        }
    }
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> Data(Int/ Torus) conversion test passed! <<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}