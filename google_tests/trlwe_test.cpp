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
    param.torusBase = 1 << 3;

    TrlweKey trlweKey {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};

    trlweKeyGen(trlweKey);

    // data gen
    DoublePolynomial plain {param.N};
    std::vector<Torus> in1(param.N);
    for (auto i = 0; i < in1.size(); i++) {
        plain.coeffs[i] = (double) genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1) / param.torusBase;
        in1[i] = doubleToTorus32(plain.coeffs[i]);
    }

    // enc
    symEncTrlweMultiSample(trlwe, trlweDft, trlweKey, in1, param.lweStdDev);
    printTrlweAB(trlwe, "trlwe");

    // decomp trlwe
    DecomposedTrlwe lhs {param};
//    for (auto l = 0; l < d.l; l++) {
//        for (auto i = 0; i < in1.size(); i++) {
//            auto decomposedMu = in1[i] << (param.torusBits - (i + 1) * param.radixBits);
//
//        }
//    }
    gadgetDecomposeTrlwe(lhs, trlwe, param);
    printDecomposedTrlweAB(lhs, "lhs");

    // first decomp
    // todo: decomposition overflow support
    Integer y = 7987;
//    vector<Integer> rhs(lhs.l);
//    decomposeOverB(rhs, y, param);
//    printArray(rhs,"rhs");
//
//    // second decomp
//    vector<DecomposedData> mid(lhs.l, DecomposedData(lhs.l));
//    for (auto i = 0; i < mid.size(); i++) {
//        gadgetDecompose(mid[i], rhs[i], param);
//        printArray(mid[i].value, "mid level " + to_string(i));
//    }

    // first recomp
    DecomposedTrlwe recomp1 {param};
//    for (auto j = 0; j < param.N; j++) {
//        for (auto r = 0; r < param.k; r++) {
//            for (auto l1 = 0; l1 < lhs.l; l1++) {
//                for (auto l2 = 0; l2 < lhs.l; l2++) {
//                    recomp1.rlwes[l1].a[r].coeffs[j] += lhs.rlwes[l1].a[r].coeffs[j] * mid[l1].value[l2];
//                }
//            }
//        }
//        for (auto l1 = 0; l1 < lhs.l; l1++) {
//            for (auto l2 = 0; l2 < lhs.l; l2++) {
//                recomp1.rlwes[l1].b.coeffs[j] += lhs.rlwes[l1].b.coeffs[j] * mid[l1].value[l2];
//            }
//        }
//    }
    for (auto j = 0; j < param.N; j++) {
        for (auto r = 0; r < param.k; r++) {
            for (auto l1 = 0; l1 < lhs.l; l1++) {
                recomp1.rlwes[l1].a[r].coeffs[j] = lhs.rlwes[l1].a[r].coeffs[j] * y;
            }
        }
        for (auto l1 = 0; l1 < lhs.l; l1++) {
            recomp1.rlwes[l1].b.coeffs[j] = lhs.rlwes[l1].b.coeffs[j] * y;
        }
    }

    printDecomposedTrlweAB(recomp1, "recomp1");

    // second recomp
    Trlwe recomp2 {param.k, param.N};
    recomposeTrlwe(recomp2, recomp1, param);
    printTrlweAB(recomp2, "recomp2");

    // dec
    DoublePolynomial dp {param.N};
    applyNttForAB(trlweDft, recomp2);
    symDecTrlwe(dp, trlweDft, trlweKey, param.torusBase);
    printArray(plain.coeffs, "p0");
    printArray(dp.coeffs, "dp");
    DoublePolynomial tv {param.N};
    for (auto i = 0; i < plain.coeffs.size(); i++) {
        tv.coeffs[i] = torus32ToDouble(doubleToTorus32(plain.coeffs[i] * y));
    }
    printArray(tv.coeffs, "p1");

    printBanner("TrlweMultConstant");
}