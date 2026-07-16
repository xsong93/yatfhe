//
// Created by Xintong Song on 2024/2/21.
//

#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(POLY, MULT0) {
    YatfheParameters param {};
    param.N = 8;
    initYatfhe(param);
    IntPolynomial polyI32A{param.N};
    IntPolynomial polyI32B001{param.N};
    IntPolynomial polyI32B010{param.N};
    IntPolynomial polyI32B011{param.N};
    IntPolynomial polyI32B100{param.N};
    IntPolynomial polyI32B101{param.N};
    IntPolynomial polyI32B110{param.N};
    IntPolynomial polyI32B111{param.N};
    IntPolynomial polyI32C{param.N};
    polyI32A.coeffs = {0,0,0,0,0,1,0,0}; // 5
    polyI32B001.coeffs = {0,1,1,1,0,0,0,0}; // =<4. i==0, i>=8-4=4
    // polyI32B010.coeffs = {0,1,0,0,0,0,0,0};
    // polyI32B011.coeffs = {0,1,1,0,0,0,0,0};
    // polyI32B100.coeffs = {1,0,0,0,0,0,0,0};
    // polyI32B101.coeffs = {1,0,1,0,0,0,0,0};
    // polyI32B110.coeffs = {1,1,0,0,0,0,0,0};
    // polyI32B111.coeffs = {1,1,1,0,0,0,0,0};

    printArray(polyI32A.coeffs, "polyI32A");
    // printArray(polyI32B001.coeffs, "polyI32B");

    multIntPolynomialModQ(polyI32C, polyI32A, polyI32B001, INT_MAX_VALUE);
    printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B010, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B011, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B100, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B101, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B110, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
    // multIntPolynomialModQ(polyI32C, polyI32A, polyI32B111, INT_MAX_VALUE);
    // printArray(polyI32C.coeffs, "polyI32C");
}

TEST(POLY, MULT) {
    YatfheParameters param {};
    param.N = 64;
    initYatfhe(param);
    IntPolynomial polyI32A{param.N};
    IntPolynomial polyI32B{param.N};
    IntPolynomial polyI32C{param.N};
    TorusPolynomial polyT32A{param.N};
    TorusPolynomial polyT32B{param.N};
    TorusPolynomial polyT32C{param.N};
    Int8Polynomial polyI8A{param.N};
    Int8Polynomial polyI8B{param.N};
    Int8Polynomial polyI8C{param.N};

    for (size_t i = 0; i < param.N; i++) {
        polyI32A.coeffs[i] = longModP(i, 4);
        polyI32B.coeffs[i] = longModP(i, 4);
        polyT32A.coeffs[i] = longModP(i, 4);
        polyT32B.coeffs[i] = longModP(i, 4);
        polyI8A.coeffs[i] = longModP(i, 4);
        polyI8B.coeffs[i] = longModP(i, 4);
    }
    printArray(polyI32A.coeffs, "polyI32A");
    printArray(polyI32B.coeffs, "polyI32B");
    printArray(polyT32A.coeffs, "polyT32A");
    printArray(polyT32B.coeffs, "polyT32B");
    printArray(polyI8A.coeffs, "polyI8A");
    printArray(polyI8B.coeffs, "polyI8B");

    multIntPolynomialAcc(polyI32C, polyI32A, polyI32B);
    multTorusPolynomialAcc(polyT32C, polyT32A, polyT32B);
    multInt8PolynomialAcc(polyI8C, polyI8A, polyI8B, 256);

    printArray(polyI32C.coeffs, "polyI32C");
    printArray(polyT32C.coeffs, "polyT32C");
    printArray(polyI8C.coeffs, "polyI8C");
    printBanner("POLY.MULT");
}

TEST(POLY, EXTERNAL_SUMPROP) {
    YatfheParameters param {};
    param.N = 64;
    initYatfhe(param);
    TorusPolynomial polyT32A{param.N};
    TorusPolynomial polyT32B{param.N};
    TorusPolynomial polyT32B1{param.N};
    TorusPolynomial polyT32C{param.N};
    TorusPolynomial polyT32C1{param.N};

    for (size_t i = 0; i < param.N; i++) {
        polyT32A.coeffs[i] = -1;
        polyT32B.coeffs[i] = genIntUniformDist(0, 1);
        polyT32B1.coeffs[i] = genIntUniformDist(0, 1);
    }
    printArray(polyT32A.coeffs, "polyT32A");
    printArray(polyT32B.coeffs, "polyT32B");
    printArray(polyT32B1.coeffs, "polyT32B1");

    TorusPolynomial BplusB1{param.N};
    TorusPolynomial addBefore{param.N};
    addTorusPolynomial(BplusB1, polyT32B, polyT32B1);
    multTorusPolynomial(addBefore, polyT32A, BplusB1);

    TorusPolynomial addAfter{param.N};
    multTorusPolynomial(polyT32C, polyT32A, polyT32B);
    multTorusPolynomial(polyT32C1, polyT32A, polyT32B1);
    addTorusPolynomial(addAfter, polyT32C, polyT32C1);

    printArray(addBefore.coeffs, "addBeforeConv");
    printArray(addAfter.coeffs, "addAfterConv");

    ASSERT_EQ(addBefore.coeffs, addAfter.coeffs);

    printBanner("POLY.EXTERNAL_SUMPROP");
}

TEST(POLY, EXTERNAL_SUMPROP2) {
    YatfheParameters param {};
    param.N = 64;
    initYatfhe(param);
    TorusPolynomial polyT32A{param.N};
    TorusPolynomial polyT32B{param.N};
    TorusPolynomial polyT32B1{param.N};
    TorusPolynomial polyT32B2{param.N};
    TorusPolynomial polyT32C{param.N};
    TorusPolynomial polyT32C1{param.N};
    TorusPolynomial polyT32C2{param.N};

    for (size_t i = 0; i < param.N; i++) {
        polyT32A.coeffs[i] = -1;
        polyT32B.coeffs[i] = genIntUniformDist(0, 1);
        polyT32B1.coeffs[i] = genIntUniformDist(0, 1);
        polyT32B2.coeffs[i] = genIntUniformDist(0, 1);
        polyT32C.coeffs[i] = genIntUniformDist(0, 1);
        polyT32C1.coeffs[i] = genIntUniformDist(0, 1);
        polyT32C2.coeffs[i] = genIntUniformDist(0, 1);
    }

    // A * (B * C + B1 * C1 + B2 * C2)
    TorusPolynomial bXc{param.N};
    multTorusPolynomialAcc(bXc, polyT32B, polyT32C);
    multTorusPolynomialAcc(bXc, polyT32B1, polyT32C1);
    multTorusPolynomialAcc(bXc, polyT32B2, polyT32C2);
    TorusPolynomial bXcXA{param.N};
    multTorusPolynomial(bXcXA, polyT32A, bXc);

    // (A * B) * C + (A * B1) * C1 + (A * B2) * C2
    TorusPolynomial aXb{param.N};
    TorusPolynomial aXb1{param.N};
    TorusPolynomial aXb2{param.N};
    TorusPolynomial aXbXc{param.N};
    TorusPolynomial aXb1Xc1{param.N};
    TorusPolynomial aXb2Xc2{param.N};
    TorusPolynomial aDbDc{param.N};
    multTorusPolynomial(aXb, polyT32A, polyT32B);
    multTorusPolynomial(aXb1, polyT32A, polyT32B1);
    multTorusPolynomial(aXb2, polyT32A, polyT32B2);
    multTorusPolynomial(aXbXc, aXb, polyT32C);
    multTorusPolynomial(aXb1Xc1, aXb1, polyT32C1);
    multTorusPolynomial(aXb2Xc2, aXb2, polyT32C2);
    addTorusPolynomial(aDbDc, aXbXc, aXb1Xc1);
    addTorusPolynomial(aDbDc, aDbDc, aXb2Xc2);

    printArray(bXcXA.coeffs, "dotBefore");
    printArray(aDbDc.coeffs, "dotAfter");

    ASSERT_EQ(bXcXA.coeffs, aDbDc.coeffs);

    printBanner("POLY.EXTERNAL_SUMPROP2");
}