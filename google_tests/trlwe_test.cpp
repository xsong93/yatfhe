//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/trglev.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

TEST(TrlweTest, TrlweEncDecSingleSampleTest) {
    YatfheParameters param {};
    param.k = 5;
    param.N = 1024;
    yatfheInit(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    Trlwe intt {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    double plain =  -1.0 / param.torusBase;
    Torus mu = doubleToTorus32(plain);

    DoublePolynomial output {param.N};
    symEncTrlweSingleSampleNtt(trlwe, trlweDft, trlweKey, mu);
//    symEncTrlweSingleSample(trlwe, trlweKey, mu);
//    printTrlweAB(trlwe, "trlwe");
//    applyInttForAB(intt, trlweDft);
//    printTrlweAB(intt, "intt");
    symDecTrlweNtt(output, trlweDft, trlweKey, param.torusBase);
//    symDecTrlweToDouble(output, trlwe, trlweKey, param.torusBase);

    cout << "mu:" << plain <<endl;
    printArray(output.coeffs, "output");
    for (auto coeff: output.coeffs) {
        ASSERT_EQ(plain, coeff);
    }
    printBanner("TrlweEncDecSingleSampleTest");
}

TEST(TrlweTest, TrlweEncDecMultiSampleTest) {
    YatfheParameters param{};
    yatfheInit(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    std::vector<int> plain(param.N);
    std::vector<Torus> in(param.N);
    for (auto i = 0; i < in.size(); i++) {
//        plain[i] = (double) genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1) / param.torusBase;
        plain[i] = genIntUniformDist(-4, 3);
        in[i] = modSwitchToTorus32(plain[i], param.torusBase);
    }
    printArray(plain, "plain");

    IntPolynomial output {param.N};
    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, in);
    printTrlweAB(trlwe, "trlwe");
    symDecTrlweToIntNtt(output, trlweDft, trlweKey, param.torusBase);
//    symEncTrlweMultiSample(trlwe, trlweKey, in);
//    symDecTrlweToDouble(output, trlwe, trlweKey, param.torusBase);

    printArray(output.coeffs, "output");
    for (auto i = 0; i < plain.size(); i++) {
        ASSERT_EQ(plain[i], output.coeffs[i]);
    }
    printBanner("TrlweEncDecMultiSampleTest");
}

TEST(TrlweTest, TRLWE_CRT_COMPOSITION) {
    YatfheParameters param{};
    yatfheInit(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    // data gen
    std::vector<int> plain(param.N);
    std::vector<Torus> in(param.N);
    for (auto i = 0; i < in.size(); i++) {
//        plain[i] = (double) genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1) / param.torusBase;
        plain[i] = genIntUniformDist(-4, 3);
        in[i] = modSwitchToTorus32(plain[i], param.torusBase);
    }
    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, in);
    printTrlweAB(trlwe, "trlwe");

    // RD
    std::vector<Trlwe8> trlweDecomp(param.d, Trlwe8{param.k, param.N});
    Trlwe trlweRecomp {param.k, param.N};
    COUNT_TIME("trlweCRTDecomp", trlweCRTDecomp(trlweDecomp, trlwe, param);)
    COUNT_TIME("trlweCRTRecomp", trlweCRTRecomp(trlweRecomp, trlweDecomp, param);)
    printTrlweAB(trlweRecomp, "trlweRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlweRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlweRecomp.b.coeffs, trlwe.b.coeffs);

    // RD 8d ver.
    Trlwe8D trlwe8D {param.k, param.N, param.d};
    Trlwe trlwe8DRecomp {param.k, param.N};
    COUNT_TIME("trlweCRTDecompNO", trlweCRTDecompNO(trlwe8D, trlwe, param);)
    COUNT_TIME("trlweCRTRecompNO", trlweCRTRecompNO(trlwe8DRecomp, trlwe8D, param);)
    printTrlweAB(trlwe8DRecomp, "trlwe8DRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlwe8DRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlwe8DRecomp.b.coeffs, trlwe.b.coeffs);

    printBanner("TRLWE_CRT_COMPOSITION");
}

TEST(TrlweTest, TRLWE_MCRT_COMPOSITION) {
    YatfheParameters param{};
    yatfheInit(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    // data gen
    std::vector<int> plain(param.N);
    std::vector<Torus> in(param.N);
    for (auto i = 0; i < in.size(); i++) {
        plain[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        in[i] = modSwitchToTorus32(plain[i], param.torusBase);
    }
    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, in);
    printTrlweAB(trlwe, "trlwe");

    // RD
    std::vector<Trlwe8> trlweDecomp(param.d, Trlwe8{param.k, param.N});
    Trlwe trlweRecomp {param.k, param.N};
    COUNT_TIME("trlweMCRTDecomp", trlweMCRTDecomp(trlweDecomp, trlwe, param);)
    COUNT_TIME("trlweMCRTToCRT", trlweMCRTToCRT(trlweDecomp, param);)
    COUNT_TIME("trlweCRTRecomp", trlweCRTRecomp(trlweRecomp, trlweDecomp, param);)
    printTrlweAB(trlweRecomp, "trlweRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlweRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlweRecomp.b.coeffs, trlwe.b.coeffs);

    printBanner("TRLWE_MCRT_COMPOSITION");
}

TEST(TrlweTest, TrlweAddSubMultiSampleTest) {
    YatfheParameters param {};
//    param.torusBase = 1 << 28;
    yatfheInit(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe1 {param.k, param.N};
    Trlwe trlwe2 {param.k, param.N};
    Trlwe trlwe3 {param.k, param.N};
    TrlweDft trlweDft1 {param.k, param.N};
    TrlweDft trlweDft2 {param.k, param.N};
    TrlweDft trlweDft3 {param.k, param.N};

    trlweKeyGen(trlweKey);

    std::vector<Integer> plain1(param.N);
    std::vector<Integer> plain2(param.N);
    std::vector<Torus> in1(param.N);
    std::vector<Torus> in2(param.N);
    for (auto i = 0; i < in1.size(); i++) {
//        plain1[i] = genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1);
//        plain2[i] = genIntUniformDist(-param.torusBase / 4 + 1, param.torusBase / 4 - 1);
        plain1[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plain2[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        in1[i] = modSwitchToTorus32(plain1[i], param.torusBase);
        in2[i] = modSwitchToTorus32(plain2[i], param.torusBase);
    }
    printArray(plain1, "plain1");
    printArray(plain2, "plain2");

    symEncTrlweMultiSampleNtt(trlwe1, trlweDft1, trlweKey, in1);
    symEncTrlweMultiSampleNtt(trlwe2, trlweDft2, trlweKey, in2);
//    symEncTrlweMultiSample(trlwe1, trlweKey, in1);
//    symEncTrlweMultiSample(trlwe2, trlweKey, in2);

    IntPolynomial output {param.N};
    TorusPolynomial torusOutput(param.N);
    IntPolynomial plainOutput(param.N);

    vectorAdd(torusOutput.coeffs, in1, in2);
    torusPolyToIntPoly(plainOutput, torusOutput, param.torusBase);
    printArray(plainOutput.coeffs, "plainOutput Add");

    trlweAddNtt(trlweDft3, trlweDft1, trlweDft2);
    symDecTrlweToIntNtt(output, trlweDft3, trlweKey, param.torusBase);
//    trlweAdd(trlwe3, trlwe1, trlwe2);
//    symDecTrlweToDouble(output, trlwe3, trlweKey, param.torusBase);
    printArray(output.coeffs, "trlweOutput Add");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_EQ(plainOutput.coeffs[i], output.coeffs[i]);
    }

    vectorSub(torusOutput.coeffs, in1, in2);
    torusPolyToIntPoly(plainOutput, torusOutput, param.torusBase);
    printArray(plainOutput.coeffs, "plainOutput Sub");

    trlweSubNtt(trlweDft3, trlweDft1, trlweDft2);
    symDecTrlweToIntNtt(output, trlweDft3, trlweKey, param.torusBase);
//    trlweSub(trlwe3, trlwe1, trlwe2);
//    symDecTrlweToDouble(output, trlwe3, trlweKey, param.torusBase);
    printArray(output.coeffs, "trlweOutput Sub");

    for (auto i = 0; i < plainOutput.N; i++) {
        ASSERT_EQ(plainOutput.coeffs[i], output.coeffs[i]);
    }

    printBanner("TrlweAddSubMultiSampleTest");
}

TEST(TrlweTest, TrlweMultLargeConstant) {
    YatfheParameters param {};
//    param.torusBase = 1 << 3;
    param.l = 4;
    yatfheInit(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
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
    Trglev trglev {param};
    trglevEncMultiSample(trglev, trlweKey, plainT, param);

    Integer y = genIntUniformDist(INT32_MIN, INT32_MAX);

    // recomp
    Trlwe recomp {param.k, param.N};
    trglevMultConst(recomp, trglev, y, param);

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
    printArray(rounded.coeffs, "rd"); // rd = p0

    printArray(plain.coeffs, "plain");
    printArray(resP.coeffs, "p1");

    for (auto i = 0 ; i < res.N; i++) {
        ASSERT_EQ(intModP(plain.coeffs[i] * y, param.torusBase), resP.coeffs[i]);
    }

    printBanner("TrlweMultLargeConstant");
}

TEST(TrlweTest, TrlweMultLargeConstantMultiLvl) {
    YatfheParameters param {};
    param.torusBase = 1 << 3;
    param.l = 4;
    param.l2 = 4;
    yatfheInit(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
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
    Trglev trglev {param};
    trglevEncMultiSample(trglev, trlweKey, plainT, param);

    Integer y = genIntUniformDist(INT32_MIN, INT32_MAX);

    Trlwe recomp2 {param.k, param.N};
    decomposedTglevMultConst(recomp2, trglev, y, param);

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

TEST(TrlweTest, TrlweDotMultLargeConstants) {
    YatfheParameters param {};
//    param.torusBase = 1 << 3;
    param.l = 4;
    param.n = 1024;
    param.N = param.n;
    yatfheInit(param);


    // prepare data
    // tlwe key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    //tlwe enc
    int pt = 1;
    Tlwe ct {param.n};
    symEncTlweSample(ct, pt, tlweKey);

    // rescale
    ScaledTlwe scaledCt {param.N * 2, param.n};
    rescaleTlweFromTorus32(scaledCt, ct);

    // todo: debug
    for (auto i = 0; i < scaledCt.n; i++) {
        scaledCt.a[i] = 2 + i;
    }


    // trlwe mult
    // trlwe key gen
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    trlweKeyGen(trlweKey);

    // data gen
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = tlweKey.s[i];
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Trglev trglev {param};
    trglevEncMultiSample(trglev, trlweKey, plainT, param);

//    Integer y = genIntUniformDist(IntMin, IntMax);

    // recomp
    Trlwe recomp {param.k, param.N};
    trglevDotMultConst(recomp, trglev, scaledCt.a, param);


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
    vectorDotMultConst(plainT.coeffs, plainT.coeffs, scaledCt.a);
    printArray(plainT.coeffs, "p0");
    printArray(res.coeffs, "re");
    printArray(rounded.coeffs, "rd"); // rd = p0

    printArray(plain.coeffs, "plain");
    printArray(resP.coeffs, "p1");

    for (auto i = 0 ; i < param.n; i++) {
        ASSERT_EQ(intModP(plain.coeffs[i] * scaledCt.a[i], param.torusBase), resP.coeffs[i]);
    }

    printBanner("TrlweMultLargeConstant");
}