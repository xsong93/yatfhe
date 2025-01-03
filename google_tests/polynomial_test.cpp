//
// Created by Xintong Song on 2024/2/21.
//

#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/tlwe.h"
#include "yautil/tool.h"

TEST(PolynomialTest, PolynomialRounding) {
    YatfheParameters param {};
//    param.torusBase = 8;
//    param.N = 1024;
//    param.n = 4;
    TorusPolynomial v {param.N};
    std::vector<Integer> d(param.N);
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    printArray(v.coeffs, "v");
    for (int i = 0; i < v.N; i++) {
        d[i] = modSwitchFromTorus32(v.coeffs[i], param.torusBase);
    }
    printArray(d, "p");

    Integer plain = 3;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TlweKey tlweKey {param.n, 0};
    lweKeyGen(tlweKey);

    Tlwe ct {param.n};
    symEncTlweSample(ct, mu, tlweKey);

    ScaledTlwe scaledCt {2 * param.N, param.n};
    rescaleTlweFromTorus32(scaledCt, ct);
    printTlweAB(ct, "ct");
    printTlweAB(scaledCt, "scaledCt");

    TorusPolynomial rpT {param.N};
    int rot = scaledCt.b;
    for (auto i = 0; i < scaledCt.n; i++) {
        rot = (rot - scaledCt.a[i] * tlweKey.s[i]) % (2 * param.N);
    }
    cout << "u*: " << rot << endl;
    torusPolynomialRotate(rpT, -rot, v);
    IntPolynomial res {param.N};
    torusPolyToIntPoly(res, rpT, param.torusBase);
    cout << "p: " << intModP(plain, param.torusBase) << endl;
    printArray(res.coeffs, "res");
    printBanner("PolynomialRounding");
}

TEST(PolynomialTest, POLY_MULT) {
    YatfheParameters param {};
    param.N = 64;
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

    polynomialMulAccNaiveI32(polyI32C, polyI32A, polyI32B);
    polynomialMulAccNaiveT32(polyT32C, polyT32A, polyT32B);
    polynomialMulAccNaiveI8(polyI8C, polyI8A, polyI8B, 256);

    printArray(polyI32C.coeffs, "polyI32C");
    printArray(polyT32C.coeffs, "polyT32C");
    printArray(polyI8C.coeffs, "polyI8C");
    printBanner("POLY_MULT");
}

TEST(PolynomialTest, POLY_EXTERNAL_SUMPROP) {
    YatfheParameters param {};
    param.N = 64;
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
    polynomialAddT32(BplusB1, polyT32B, polyT32B1);
    polynomialMulNaiveT32(addBefore, polyT32A, BplusB1);

    TorusPolynomial addAfter{param.N};
    polynomialMulNaiveT32(polyT32C, polyT32A, polyT32B);
    polynomialMulNaiveT32(polyT32C1, polyT32A, polyT32B1);
    polynomialAddT32(addAfter, polyT32C, polyT32C1);

    printArray(addBefore.coeffs, "addBeforeConv");
    printArray(addAfter.coeffs, "addAfterConv");

    ASSERT_EQ(addBefore.coeffs, addAfter.coeffs);

    printBanner("POLY_EXTERNAL_SUMPROP");
}

TEST(PolynomialTest, POLY_EXTERNAL_SUMPROP2) {
    YatfheParameters param {};
    param.N = 64;
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
    polynomialMulAccNaiveT32(bXc, polyT32B, polyT32C);
    polynomialMulAccNaiveT32(bXc, polyT32B1, polyT32C1);
    polynomialMulAccNaiveT32(bXc, polyT32B2, polyT32C2);
    TorusPolynomial bXcXA{param.N};
    polynomialMulNaiveT32(bXcXA, polyT32A, bXc);

    // (A * B) * C + (A * B1) * C1 + (A * B2) * C2
    TorusPolynomial aXb{param.N};
    TorusPolynomial aXb1{param.N};
    TorusPolynomial aXb2{param.N};
    TorusPolynomial aXbXc{param.N};
    TorusPolynomial aXb1Xc1{param.N};
    TorusPolynomial aXb2Xc2{param.N};
    TorusPolynomial aDbDc{param.N};
    polynomialMulNaiveT32(aXb, polyT32A, polyT32B);
    polynomialMulNaiveT32(aXb1, polyT32A, polyT32B1);
    polynomialMulNaiveT32(aXb2, polyT32A, polyT32B2);
    polynomialMulNaiveT32(aXbXc, aXb, polyT32C);
    polynomialMulNaiveT32(aXb1Xc1, aXb1, polyT32C1);
    polynomialMulNaiveT32(aXb2Xc2, aXb2, polyT32C2);
    polynomialAddT32(aDbDc, aXbXc, aXb1Xc1);
    polynomialAddT32(aDbDc, aDbDc, aXb2Xc2);

    printArray(bXcXA.coeffs, "dotBefore");
    printArray(aDbDc.coeffs, "dotAfter");

    ASSERT_EQ(bXcXA.coeffs, aDbDc.coeffs);

    printBanner("POLY_EXTERNAL_SUMPROP2");
}