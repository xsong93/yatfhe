//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
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
    symEncTrlweSingleSample(trlwe, trlweDft, trlweKey, mu, param.lweStdDev);
    printTrlweAB(trlwe, "trlwe");
    applyInttForAB(intt, trlweDft);
    printTrlweAB(intt, "intt");
    symDecTrlwe(output, trlweDft, trlweKey, param.torusBase);

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
    symEncTrlweMultiSample(trlwe, trlweDft, trlweKey, in, param.lweStdDev);
    symDecTrlwe(output, trlweDft, trlweKey, param.torusBase);

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

    symEncTrlweMultiSample(trlwe1, trlweDft1, trlweKey, in1, param.lweStdDev);
    symEncTrlweMultiSample(trlwe2, trlweDft2, trlweKey, in2, param.lweStdDev);

    DoublePolynomial output {param.N};
    TorusPolynomial torusOutput(param.N);
    DoublePolynomial  plainOutput(param.N);

    vectorAdd(torusOutput.coeffs, in1, in2);
    torusPolyToDoublePoly(plainOutput, torusOutput);
    printArray(plainOutput.coeffs, "plainOutput Add");

    trlweAddNtt(trlweDft3, trlweDft1, trlweDft2);
    symDecTrlwe(output, trlweDft3, trlweKey, param.torusBase);
    printArray(output.coeffs, "output Add");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_NEAR(plainOutput.coeffs[i], output.coeffs[i], 0.001);
    }

    vectorSub(torusOutput.coeffs, in1, in2);
    torusPolyToDoublePoly(plainOutput, torusOutput);
    printArray(plainOutput.coeffs, "plainOutput Sub");

    trlweSubNtt(trlweDft3, trlweDft1, trlweDft2);
    symDecTrlwe(output, trlweDft3, trlweKey, param.torusBase);
    printArray(output.coeffs, "output Sub");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_NEAR(plainOutput.coeffs[i], output.coeffs[i], 0.001);
    }

    printBanner("TrlweAddSubMultiSampleTest");
}