//
// Created by Xintong Song on 2024/2/29.
//
#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yautil/tool.h"

TEST(intTorusTest, intTorusTest) {
    int N = 100;
    IntPolynomial intPoly(N);
    TorusPolynomial tPoly(N);
    int i = -50;
    for (auto& item : intPoly.coeffs) {
        item = i++;
    }
    printArray(intPoly.coeffs, "intPoly");
    for (auto j = 0; j < N; j++) {
        tPoly.coeffs[j] = modSwitchToTorus32(intPoly.coeffs[j], 8);
    }
    printArray(tPoly.coeffs, "tPoly");
    for (auto j = 0; j < N; j++) {
        tPoly.coeffs[j] = modSwitchToTorus32(intPoly.coeffs[j], 8);
    }
}