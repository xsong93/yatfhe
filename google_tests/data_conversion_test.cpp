//
// Created by Xintong Song on 2024/2/29.
//
#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yautil/tool.h"


TEST(doubleTorusTest, doubleTorusTest) {
    int N = 100;
    vector<double> doublePoly(N);
    TorusPolynomial tPoly(N);
    vector<double> resPoly(N);
    double i = -0.5;
    for (auto& item : doublePoly) {
        item = i;
        i += 0.01;
    }
    printArray(doublePoly, "doublePoly");
    for (auto j = 0; j < N; j++) {
        tPoly.coeffs[j] = doubleToTorus32(doublePoly[j]);
    }
    printArray(tPoly.coeffs, "torusPoly");
    for (auto j = 0; j < N; j++) {
        resPoly[j] = torus32ToDouble(tPoly.coeffs[j]);
    }
    printArray(resPoly, "resPoly");
    for (auto j = 0; j < N; j++) {
        EXPECT_NEAR(resPoly[j], doublePoly[j], 1e-3);
    }
    std::cout << ">>>>>>>> Data(double/ Torus) conversion test passed! <<<<<<<<" << std::endl;
}

TEST(intTorusTest, intTorusTest) {
    int N = 100;
    IntPolynomial intPoly(N);
    TorusPolynomial tPoly(N);
    IntPolynomial resPoly(N);
    int i = -50;
    for (auto& item : intPoly.coeffs) {
        item = i++;
    }
    printArray(intPoly.coeffs, "intPoly");
    for (auto j = 0; j < N; j++) {
        tPoly.coeffs[j] = modSwitchToTorus32(intPoly.coeffs[j], N);
    }
    printArray(tPoly.coeffs, "torusPoly");
    for (auto j = 0; j < N; j++) {
        resPoly.coeffs[j] = modSwitchFromTorus32(tPoly.coeffs[j], N);
    }
    printArray(resPoly.coeffs, "resPoly");
    for (auto j = 0; j < N; j++) {
        EXPECT_EQ(resPoly.coeffs[j], intPoly.coeffs[j]);
    }
    std::cout << ">>>>>>>> Data(Int/ Torus) conversion test passed! <<<<<<<<" << std::endl;
}