//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/tglev.h"
#include "yautil/tool.h"

TEST(TrlweEncDecSingleSampleTest, TrlweEncDecSingleSampleTest) {
    YatfheParameters param {};
    param.k = 5;
    param.N = 1024;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    Trlwe intt {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    double plain =  -1.0 / param.torusBase;
    Torus mu = doubleToTorus32(plain);

    DoublePolynomial output {param.N};
//    symEncTrlweSingleSampleNtt(trlwe, trlweDft, trlweKey, mu, param.rlweStdDev);
    symEncTrlweSingleSample(trlwe, trlweKey, mu, param.rlweStdDev);
//    printTrlweAB(trlwe, "trlwe");
//    applyInttForAB(intt, trlweDft);
//    printTrlweAB(intt, "intt");
//    symDecTrlweNtt(output, trlweDft, trlweKey, param.torusBase);
    symDecTrlwe(output, trlwe, trlweKey, param.torusBase);

    cout << "mu:" << plain <<endl;
    printArray(output.coeffs, "output");
    for (auto coeff: output.coeffs) {
        ASSERT_EQ(plain, coeff);
    }
    printBanner("TrlweEncDecSingleSampleTest");
}

TEST(TrlweEncDecMultiSampleTest, TrlweEncDecMultiSampleTest) {
    const YatfheParameters param {};

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    std::vector<double> plain(param.N);
    std::vector<Torus> in(param.N);
    for (auto i = 0; i < in.size(); i++) {
        plain[i] = (double) genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1) / param.torusBase;
        in[i] = doubleToTorus32(plain[i]);
    }
    printArray(plain, "plain");

    DoublePolynomial output {param.N};
//    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, in, param.rlweStdDev);
//    symDecTrlweNtt(output, trlweDft, trlweKey, param.torusBase);
    symEncTrlweMultiSample(trlwe, trlweKey, in, param.rlweStdDev);
    symDecTrlwe(output, trlwe, trlweKey, param.torusBase);

    printArray(output.coeffs, "output");
    for (auto i = 0; i < plain.size(); i++) {
        ASSERT_EQ(plain[i], output.coeffs[i]);
    }
    printBanner("TrlweEncDecMultiSampleTest");
}

TEST(TrlweAddSubMultiSampleTest, TrlweAddSubMultiSampleTest) {
    YatfheParameters param {};
    param.torusBase = 1 << 28;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe1 {param.k, param.N};
    Trlwe trlwe2 {param.k, param.N};
    Trlwe trlwe3 {param.k, param.N};
    TrlweDft trlweDft1 {param.k, param.N};
    TrlweDft trlweDft2 {param.k, param.N};
    TrlweDft trlweDft3 {param.k, param.N};

    trlweKeyGen(trlweKey);

    std::vector<double> plain1(param.N);
    std::vector<double> plain2(param.N);
    std::vector<Torus> in1(param.N);
    std::vector<Torus> in2(param.N);
    for (auto i = 0; i < in1.size(); i++) {
        plain1[i] = (double) genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1) / param.torusBase;
        plain2[i] = (double) genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1) / param.torusBase;
        in1[i] = doubleToTorus32(plain1[i]);
        in2[i] = doubleToTorus32(plain2[i]);
    }
    printArray(plain1, "plain1");
    printArray(plain2, "plain2");

//    symEncTrlweMultiSampleNtt(trlwe1, trlweDft1, trlweKey, in1, param.rlweStdDev);
//    symEncTrlweMultiSampleNtt(trlwe2, trlweDft2, trlweKey, in2, param.rlweStdDev);
    symEncTrlweMultiSample(trlwe1, trlweKey, in1, param.rlweStdDev);
    symEncTrlweMultiSample(trlwe2, trlweKey, in2, param.rlweStdDev);

    DoublePolynomial output {param.N};
    TorusPolynomial torusOutput(param.N);
    DoublePolynomial plainOutput(param.N);

    vectorAdd(torusOutput.coeffs, in1, in2);
    torusPolyToDoublePoly(plainOutput, torusOutput);
    printArray(plainOutput.coeffs, "plainOutput Add");

//    trlweAddNtt(trlweDft3, trlweDft1, trlweDft2);
//    symDecTrlweNtt(output, trlweDft3, trlweKey, param.torusBase);
    trlweAdd(trlwe3, trlwe1, trlwe2);
    symDecTrlwe(output, trlwe3, trlweKey, param.torusBase);
    printArray(output.coeffs, "output Add");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_NEAR(plainOutput.coeffs[i], output.coeffs[i], 0.001);
    }

    vectorSub(torusOutput.coeffs, in1, in2);
    torusPolyToDoublePoly(plainOutput, torusOutput);
    printArray(plainOutput.coeffs, "plainOutput Sub");

//    trlweSubNtt(trlweDft3, trlweDft1, trlweDft2);
//    symDecTrlweNtt(output, trlweDft3, trlweKey, param.torusBase);
    trlweSub(trlwe3, trlwe1, trlwe2);
    symDecTrlwe(output, trlwe3, trlweKey, param.torusBase);
    printArray(output.coeffs, "output Sub");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_NEAR(plainOutput.coeffs[i], output.coeffs[i], 0.001);
    }

    printBanner("TrlweAddSubMultiSampleTest");
}

TEST(TrlweMultLargeConstant, TrlweMultLargeConstant) {
    YatfheParameters param {};
    param.torusBase = 1 << 3;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};

    trlweKeyGen(trlweKey);

    // data gen
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Tglev tglev {param};
    tglevEncMultiSample(tglev, trlweKey, plainT, param);

    Integer y = genIntUniformDist(IntMin, IntMax);

    // recomp
    Trlwe recomp {param.k, param.N};
    tglevMultConst(recomp, tglev, y, param);

    // dec
    TorusPolynomial res {param.N};
    TorusPolynomial rounded {param.N};
    IntPolynomial resP {param.N};
//    applyNttForAB(trlweDft, recomp);
//    symDecTrlweWoRoundingNtt(res, trlweDft, trlweKey);
    symDecTrlweWoRounding(res, recomp, trlweKey);


    for (auto i = 0 ; i < res.N; i++) {
        rounded.coeffs[i] = roundTorusError(res.coeffs[i], param.torusBase);
        resP.coeffs[i] = modSwitchFromTorus32(rounded.coeffs[i], param.torusBase);
    }

    printArray(plainT.coeffs, "plainT");
    vectorMultConst(plainT.coeffs, plainT.coeffs, y);
    printArray(plainT.coeffs, "p0");
    printArray(res.coeffs, "re");
    printArray(rounded.coeffs, "rd");

    printArray(plain.coeffs, "plain");
    printArray(resP.coeffs, "p1");

    for (auto i = 0 ; i < res.N; i++) {
        ASSERT_EQ(intModP(plain.coeffs[i] * y, param.torusBase), resP.coeffs[i]);
    }

    printBanner("TrlweMultLargeConstant");
}

TEST(TrlweMultLargeConstantMultiLvl, TrlweMultLargeConstantMultiLvl) {
    YatfheParameters param {};
    param.torusBase = 1 << 3;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};

    trlweKeyGen(trlweKey);

    // data gen
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Tglev tglev {param};
    tglevEncMultiSample(tglev, trlweKey, plainT, param);

    Integer y = genIntUniformDist(IntMin, IntMax);

    Trlwe recomp2 {param.k, param.N};
    decomposedTglevMultConst(recomp2, tglev, y, param);

    // dec
    TorusPolynomial res {param.N};
    TorusPolynomial rounded {param.N};
    IntPolynomial resP {param.N};
    symDecTrlweWoRounding(res, recomp2, trlweKey);

    for (auto i = 0 ; i < res.N; i++) {
        rounded.coeffs[i] = roundTorusError(res.coeffs[i], param.torusBase);
        resP.coeffs[i] = modSwitchFromTorus32(rounded.coeffs[i], param.torusBase);
    }

    printArray(plainT.coeffs, "plainT");
    vectorMultConst(plainT.coeffs, plainT.coeffs, y);
    printArray(plainT.coeffs, "p0");
    printArray(res.coeffs, "re");
    printArray(rounded.coeffs, "rd");

    printArray(plain.coeffs, "plain");
    printArray(resP.coeffs, "p1");

    for (auto i = 0 ; i < res.N; i++) {
        ASSERT_EQ(intModP(plain.coeffs[i] * y, param.torusBase), resP.coeffs[i]);
    }

    printBanner("TrlweMultLargeConstantMultiLvl");
}