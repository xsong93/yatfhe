//
// Created by Xintong Song on 2024/2/29.
//
#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(DATA_CONVERT, DOUBLE_TORUS) {
    YatfheParameters param {};
    initYatfhe(param);
    auto N = 100;

    DoublePolynomial doublePoly{N};
    TorusPolynomial tPoly{N};
    TorusPolynomial t2Poly{N};
    DoublePolynomial resPoly{N};

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
        initCoeffsViaUniformDistribution(tPoly.coeffs, TORUS_MIN, TORUS_MAX);
//        printArray(tPoly.coeffs, "torusPoly");
        torusPolyToDoublePoly(resPoly, tPoly);
//        printArray(resPoly.coeffs, "resPoly");
        doublePolyToTorusPoly(t2Poly, resPoly);
//        printArray(t2Poly.coeffs, "t2Poly");
        for (auto j = 0; j < N; j++) {
            EXPECT_NEAR(tPoly.coeffs[j], t2Poly.coeffs[j], 1e3);
        }
    }
    printBanner("DATA_CONVERT.DOUBLE_TORUS");
}

TEST(DATA_CONVERT, INT_TORUS) {
    YatfheParameters param {};
    initYatfhe(param);
    int N = 100;

    IntPolynomial intPoly{N};
    TorusPolynomial tPoly{N};
    IntPolynomial resPoly{N};
    TorusPolynomial resTPoly{N};
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
        initCoeffsViaUniformDistribution(tPoly.coeffs, TORUS_MIN, TORUS_MAX);
//        printArray(tPoly.coeffs, "torusPoly");
        torusPolyToIntPoly(resPoly, tPoly, 2 * N);
//        printArray(resPoly.coeffs, "resPoly");
        intPolyToTorusPoly(resTPoly, resPoly, 2 * N);
//        printArray(resTPoly.coeffs, "torusPoly");
        roundErrorTorusPoly(tPoly, 2 * N);
        for (auto j = 0; j < N; j++) {
            EXPECT_EQ(resTPoly.coeffs[j], tPoly.coeffs[j]);
        }
    }
    printBanner("DATA_CONVERT.INT_TORUS");
}

TEST(DATA_CONVERT, SCALE_TLWE) {
    YatfheParameters param {};
    param.torusBase = 8;

    initYatfhe(param);

    for (size_t t = 0; t < 100; t++) {
        TlweKey tlweKey{param};

        genTlweKey(tlweKey);

        int plain = 2;
        Torus mu = modSwitchToTorus32(plain, param.torusBase);

        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);

        ScaledTlwe inputModN2{param.N * 2, param.n};
        rescaleTlweToNewMod(inputModN2, input);

        // printTlweAB(input, "input");
        // printTlweAB(inputModN2, "inputModN2");

        for (size_t i = 0; i < param.n; i++) {
            ASSERT_LT(inputModN2.a[i], param.N);
            ASSERT_GE(inputModN2.a[i], -param.N);
        }
        ASSERT_LE(inputModN2.b, param.N);
        ASSERT_GE(inputModN2.b, -param.N);
    }
    printBanner("DATA_CONVERT.SCALE_TLWE");
}

TEST(DATA_CONVERT, ROUND_TORUS) {
    cout << modSwitchToTorus32(2, 8) << endl;
    cout << modSwitchToTorus32(3, 8) << endl;
    cout << modSwitchToTorus32(4, 8) << endl;

    cout << modSwitchFromTorus32(631605433, 8) << endl;

    cout << modSwitchFromTorus32(631605434, 8) << endl;
    cout << modSwitchFromTorus32(842140578, 8) << endl;
    cout << modSwitchFromTorus32(1052675722, 8) << endl;

    cout << modSwitchFromTorus32(1052675723, 8) << endl;

    cout << roundTorus32Error(631605433, 8) << endl;
    cout << roundTorus32Error(631605434, 8) << endl;
    cout << roundTorus32Error(842140578, 8) << endl;
    cout << roundTorus32Error(1052675722, 8) << endl;
    cout << roundTorus32Error(1052675723, 8) << endl;
}