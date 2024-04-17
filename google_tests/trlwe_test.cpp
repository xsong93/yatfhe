//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/gadget_decomposition.h"
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

TEST(TrlweMultConstant, TrlweMultConstant) {
    YatfheParameters param {};
    param.torusBase = 1 << 28;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe1 {param.k, param.N};
    TrlweDft trlweDft1 {param.k, param.N};

    trlweKeyGen(trlweKey);

    std::vector<double> plain1(param.N);
    std::vector<Torus> in1(param.N);
    for (auto i = 0; i < in1.size(); i++) {
        plain1[i] = (double) genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1) / param.torusBase;
        in1[i] = doubleToTorus32(plain1[i]);
    }
    printArray(plain1, "plain1");

    symEncTrlweMultiSample(trlwe1, trlweDft1, trlweKey, in1, param.lweStdDev);
    printTrlweAB(trlwe1, "t1");

    DecomposedTrlwe d {param};
//    for (auto l = 0; l < d.l; l++) {
//        for (auto i = 0; i < in1.size(); i++) {
//            auto decomposedMu = in1[i] << (param.torusBits - (i + 1) * param.radixBits);
//
//        }
//    }
    gadgetDecomposeTrlwe(d, trlwe1, param);
    printDecomposedTrlweAB(d, "d");

    // todo: fix overflow issue
    Integer y = 5;
    vector<Integer> yOr(d.l);
    decomposeOverB(yOr, y, param);
    printArray(yOr,"yOr");
    vector<DecomposedData> dy2(d.l, DecomposedData(d.l));
    for (auto i = 0; i < dy2.size(); i++) {
        gadgetDecompose(dy2[i], yOr[i], param);
        printArray(dy2[i].value, "dy2 level " + to_string(i));
    }


//    Trlwe out {param.k, param.N};
    DecomposedTrlwe out {param};
    for (auto j = 0; j < param.N; j++) {
        for (auto r = 0; r < param.k; r++) {
            for (auto l1 = 0; l1 < d.l; l1++) {
                for (auto l2 = 0; l2 < d.l; l2++) {
                    out.rlwes[l1].a[r].coeffs[j] += d.rlwes[l1].a[r].coeffs[j] * dy2[l1].value[l2];
                }
            }
        }
        for (auto l1 = 0; l1 < d.l; l1++) {
            for (auto l2 = 0; l2 < d.l; l2++) {
                out.rlwes[l1].b.coeffs[j] += d.rlwes[l1].b.coeffs[j] * dy2[l1].value[l2];
            }
        }
    }

    printDecomposedTrlweAB(out, "out");
    Trlwe fi {param.k, param.N};
    recomposeTrlwe(fi, out, param);
    printTrlweAB(fi, "fi");

    DoublePolynomial dp {param.N};
    TrlweDft trlweDft {param.k, param.N};
    applyNttForAB(trlweDft, fi);
    symDecTrlwe(dp, trlweDft, trlweKey, param.torusBase);
    printArray(dp.coeffs, "dp");
    DoublePolynomial tv {param.N};
    for (auto i = 0; i < plain1.size(); i++) {
        tv.coeffs[i] = plain1[i] * y;
    }
    printArray(tv.coeffs, "tv");

    printBanner("TrlweMultConstant");
}