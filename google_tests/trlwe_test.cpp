//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

TEST(TrlweEncDecSingleSampleTest, TrlweEncDecSingleSampleTest) {
    const YatfheParameters param {};

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    double plain =  -1.0 / param.torusBase;
    Torus mu = doubleToTorus32(plain);

    DoublePolynomial output {param.N};
    symEncTrlweSingleSample(trlwe, trlweDft, trlweKey, mu, param.lweStdDev);
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