//
// Created by Xintong Song on 2024/3/19.
//
#include <thread>

#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/trglev.h"
#include "yatfhe/crt.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

TEST(TrlweTest, TRLWE_RED_TEST) {
    YatfheParameters param {};
    initYatfhe(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    Trlwe intt {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

    double plain =  -1.0 / param.torusBase;
    Torus mu = doubleToTorus32(plain);

    DoublePolynomial output {param.N};
    symEncTrlweSingleSample(trlwe, trlweKey, mu);

    // cut off
    int thres = 16;
    for (size_t i = 0; i < param.k; i++) {
        for (size_t j = 0; j < param.N; j++) {
            auto tmp = trlwe.a[i].coeffs[j] >> thres;
            trlwe.a[i].coeffs[j] = tmp << thres;
        }
    }
    for (size_t j = 0; j < param.N; j++) {
        auto tmp = trlwe.b.coeffs[j] >> thres;
        trlwe.b.coeffs[j] = tmp << thres;
    }

    symDecTrlweToDouble(output, trlwe, trlweKey, param.torusBase);

    cout << "mu:" << plain <<endl;
    printArray(output.coeffs, "output");
    for (auto coeff: output.coeffs) {
        ASSERT_EQ(plain, coeff);
    }
    printBanner("TRLWE_RED_TEST");
}

TEST(TrlweTest, TrlweEncDecSingleSampleTest) {
    YatfheParameters param {};
    param.k = 5;
    param.N = 1024;
    initYatfhe(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    Trlwe intt {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

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
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

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

    printArray(output.coeffs, "output");
    for (auto i = 0; i < plain.size(); i++) {
        ASSERT_EQ(plain[i], output.coeffs[i]);
    }
    printBanner("TrlweEncDecMultiSampleTest");
}

TEST(TrlweTest, TRLWE_ENCS) {
    YatfheParameters param{};
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

    std::vector<IntPolynomial> plain(param.k, TorusPolynomial(param.N));
    std::vector<TorusPolynomial> in(param.k, TorusPolynomial(param.N));
    std::vector<TorusPolynomial> negativeS(param.k, TorusPolynomial(param.N));
    for (size_t i = 0; i < param.k; i++) {
        for (size_t j = 0; j < param.N; j++) {
            plain[i].coeffs[j] = 0;
            in[i].coeffs[j] = modSwitchToTorus32(plain[i].coeffs[j], param.torusBase);
            negativeS[i].coeffs[j] = -trlweKey.s[i].coeffs[j];
        }
    }
    int pos = 0;
    int val = 1;
    in[0].coeffs[pos] = modSwitchToTorus32(val, param.torusBase);
    in[1].coeffs[pos] = modSwitchToTorus32(val, param.torusBase);

    std::vector<TorusPolynomial> sXm(param.k, TorusPolynomial(param.N));
    for (size_t i = 0; i < param.k; i++) {
        multTorusPolynomial(sXm[i], in[i], negativeS[i]);
    }
    printArray(sXm[0].coeffs, "-sXm 0");
    printArray(sXm[1].coeffs, "-sXm 1");

    // RLWE(-sm) = RLWE(0) - m * (-1, 0) = (a + m, as + e) [DM'21]
    std::vector<Trlwe> encSxM(param.k, Trlwe{param.k, param.N});
    for (size_t k = 0; k < param.k; k++) {
        symEncTrlweSingleSample(encSxM[k], trlweKey, 0);
        for (size_t i = 0; i < param.k; i++) {
            if (i == k) {
                addTorusPolynomial(encSxM[k].a[i], encSxM[k].a[i], in[i]);
            }
        }
    }

    std::vector<TorusPolynomial> res(param.k, TorusPolynomial(param.N));
    symDecTrlweToTorus(res[0], encSxM[0], trlweKey, param.torusBase);
    symDecTrlweToTorus(res[1], encSxM[1], trlweKey, param.torusBase);
    printArray(res[0].coeffs, "res 0");
    printArray(res[1].coeffs, "res 1");

    ASSERT_EQ(res[0].coeffs, sXm[0].coeffs);
    ASSERT_EQ(res[1].coeffs, sXm[1].coeffs);

    printBanner("TRLWE_ENCS");
}

TEST(TrlweTest, TRLWE_ROT) {
    YatfheParameters param{};
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

    std::vector<int> plain(param.N);
    std::vector<Torus> in(param.N);
    for (auto i = 0; i < in.size(); i++) {
//        plain[i] = (double) genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1) / param.torusBase;
        plain[i] = genIntUniformDist(-4, 3);
        in[i] = modSwitchToTorus32(plain[i], param.torusBase);
    }
    printArray(plain, "plain");

    IntPolynomial output {param.N};

    int rotN = -1;
    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, in);

    Trlwe rot{param.k, param.N};
    COUNT_TIME("rotateTrlwe", rotateTrlwe(rot, trlwe, rotN);)
    symDecTrlweToInt(output, rot, trlweKey, param.torusBase);
    printArray(output.coeffs, "output");
    TorusPolynomial res{param.N};
    rotateTorusPolynomial(res, -rotN, output);
    printArray(res.coeffs, "res");
    ASSERT_EQ(plain, res.coeffs);

    TrlweDft rotDft {param.k, param.N};
    COUNT_TIME("rotateTrlweNtt", rotateTrlweNtt(rotDft, trlweDft, rotN);)
    symDecTrlweToIntNtt(output, rotDft, trlweKey, param.torusBase);
    printArray(output.coeffs, "outputDft");
    rotateTorusPolynomial(res, -rotN, output);
    printArray(res.coeffs, "resDft");
    ASSERT_EQ(plain, res.coeffs);

    printBanner("TRLWE_ROT");
}

TEST(TrlweTest, TRLWE_CRT_COMPOSITION) {
    YatfheParameters param{};
    param.q = Q_CRT;
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

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
    COUNT_TIME("trlweCRTDecomp", decompTrlweCrt(trlweDecomp, trlwe, param);)
    COUNT_TIME("trlweCRTRecomp", recompTrlweCrt(trlweRecomp, trlweDecomp, param);)
    printTrlweAB(trlweRecomp, "trlweRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlweRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlweRecomp.b.coeffs, trlwe.b.coeffs);

    // RD 8d ver.
    Trlwe8D trlwe8D {param.k, param.N, param.d};
    Trlwe trlwe8DRecomp {param.k, param.N};
    COUNT_TIME("trlweCRTDecompNO", decompTrlweCrtNO(trlwe8D, trlwe, param);)
    COUNT_TIME("trlweCRTRecompNO", recompTrlweCrtNO(trlwe8DRecomp, trlwe8D, param);)
    printTrlweAB(trlwe8DRecomp, "trlwe8DRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlwe8DRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlwe8DRecomp.b.coeffs, trlwe.b.coeffs);

    printBanner("TRLWE_CRT_COMPOSITION");
}

TEST(TrlweTest, TRLWE_MCRT_COMPOSITION) {
    YatfheParameters param{};
    param.q = Q_CRT;
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

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
    COUNT_TIME("trlweMCRTDecomp", decompTrlweMcrt(trlweDecomp, trlwe, param);)
    COUNT_TIME("trlweMCRTToCRT", trlweMcrtToCrt(trlweDecomp, param);)
    COUNT_TIME("trlweCRTRecomp", recompTrlweCrt(trlweRecomp, trlweDecomp, param);)
    printTrlweAB(trlweRecomp, "trlweRecomp");
    for (auto i = 0; i < param.k; i++) {
        ASSERT_EQ(trlweRecomp.a[i].coeffs, trlwe.a[i].coeffs);
    }
    ASSERT_EQ(trlweRecomp.b.coeffs, trlwe.b.coeffs);

    printBanner("TRLWE_MCRT_COMPOSITION");
}

TEST(TrlweTest, TRLWE_APPROX_CRT_COMPOSITION) {
    YatfheParameters param{};
    param.q = Q_CRT;
    initYatfhe(param);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};
    genTrlweKey(trlweKey);

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
    std::vector<Trlwe8> trlweGadDecomp(param.d, Trlwe8{param.k, param.N});
    Trlwe trlweRecomp {param.k, param.N};
    COUNT_TIME("trlweMCRTDecomp", decompTrlweMcrt(trlweDecomp, trlwe, param);)
    COUNT_TIME("trlweApproxCRTDecomp", decompTrlweApproxCrt(trlweGadDecomp, trlweDecomp, param);)
    COUNT_TIME("trlweApproxCRTRecomp", recompTrlweApproxCrt(trlweRecomp, trlweGadDecomp, param);)
    printTrlweAB(trlweRecomp, "trlweRecomp");
    int errA = 0;
    int errB = 0;
    for (auto i = 0; i < param.k; i++) {
        for (auto j = 0; j < param.N; j++) {
            errA = std::max(errA, std::abs(trlweRecomp.a[i].coeffs[j] - trlwe.a[i].coeffs[j]));
        }
    }
    for (auto j = 0; j < param.N; j++) {
        errB = std::max(errB, std::abs(trlweRecomp.b.coeffs[j] - trlwe.b.coeffs[j]));
    }
    int maxErr = param.dl * param.qLow / 2;
    cout <<"errA: " << errA << ", errB: "<< errB << ", maxE: " << maxErr << endl << endl;
    ASSERT_LE(errA, maxErr);
    ASSERT_LE(errB, maxErr);

    // Dec
    TorusPolynomial ori{param.N};
    TorusPolynomial gd{param.N};
    symDecTrlweToInt(ori, trlwe, trlweKey, param.torusBase);
    symDecTrlweToInt(gd, trlweRecomp, trlweKey, param.torusBase);
    printArray(ori.coeffs, "ori");
    printArray(gd.coeffs, "gd");
    ASSERT_EQ(ori.coeffs, gd.coeffs);

    printBanner("TRLWE_APPROX_CRT_COMPOSITION");
}

TEST(TrlweTest, TrlweAddSubMultiSampleTest) {
    YatfheParameters param {};
//    param.torusBase = 1 << 28;
    initYatfhe(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe1 {param.k, param.N};
    Trlwe trlwe2 {param.k, param.N};
    Trlwe trlwe3 {param.k, param.N};
    TrlweDft trlweDft1 {param.k, param.N};
    TrlweDft trlweDft2 {param.k, param.N};
    TrlweDft trlweDft3 {param.k, param.N};

    genTrlweKey(trlweKey);

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

    addTrlweNtt(trlweDft3, trlweDft1, trlweDft2);
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

    subTrlweNtt(trlweDft3, trlweDft1, trlweDft2);
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
    initYatfhe(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};

    genTrlweKey(trlweKey);

    // data gen
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Trglev trglev {param};
    encTrglevMultiSample(trglev, trlweKey, plainT, param);

    Integer y = 3;

    // recomp
    Trlwe recomp {param.k, param.N};
    multTrglevWithConst(recomp, trglev, y, param);

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

TEST(TrlweTest, TRLWE_MULT_DECOMP) {
    YatfheParameters param{};
//    param.torusBase = 1 << 3;
    param.l = 4;
    initYatfhe(param);

    TrlweKey trlweKey{param.k, param.N, param.rlweStdDev};
    genTrlweKey(trlweKey);

    // data gen
    IntPolynomial plain{param.N}; // Z/pZ
    vector<Torus> plainT(param.N);
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Trlwe trlwe{param.k, param.N};
    symEncTrlweMultiSample(trlwe, trlweKey, plainT);

    // decomp mult
    Integer a = 3;
    DecomposedData da{param.ksLevel};
    gadgetDecompose(da, a, param);
    printArray(da.value, "da");

    DecomposedTrlwe trlwes{param};
    for (size_t l = 0; l < param.l; l++) {
        for (size_t i = 0; i < param.k; i++) {
            for (size_t j = 0; j < param.N; j++) {
                trlwes.trlwes[l].a[i].coeffs[j] = multTorus(trlwe.a[i].coeffs[j], da.value[l]);
            }
        }
        for (size_t j = 0; j < param.N; j++) {
            trlwes.trlwes[l].b.coeffs[j] = multTorus(trlwe.b.coeffs[j], da.value[l]);
        }
    }
    printDecomposedTrlweAB(trlwes, "trlwes");

    // recomp
    Trlwe recomp{param.k, param.N};
    recomposeTrlwe(recomp, trlwes, param);

    // dec
    TorusPolynomial res{param.N};
    TorusPolynomial rounded {param.N};
    IntPolynomial resP{param.N};
    symDecTrlweWoRounding(res, recomp, trlweKey);


    for (auto i = 0 ; i < res.N; i++) {
        rounded.coeffs[i] = roundTorusError(res.coeffs[i], param.torusBase);
        resP.coeffs[i] = modSwitchFromTorus32(rounded.coeffs[i], param.torusBase);
    }

    printArray(plainT, "plainT");
    vectorMultConstModQ(plainT, plainT, a, TORUS_Q);
    printArray(plainT, "p0");
    printArray(res.coeffs, "re");
    printArray(rounded.coeffs, "rd"); // rd = p0

    printArray(plain.coeffs, "plain");
    printArray(resP.coeffs, "p1");

    for (auto i = 0 ; i < res.N; i++) {
        ASSERT_EQ(intModP(plain.coeffs[i] * a, param.torusBase), resP.coeffs[i]);
    }

    printBanner("TRLWE_MULT_DECOMP");
}

TEST(TrlweTest, TrlweMultLargeConstantMultiLvl) {
    YatfheParameters param {};
    param.torusBase = 1 << 3;
    param.l = 4;
    param.l2 = 4;
    initYatfhe(param);

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    Trlwe trlwe {param.k, param.N};
    TrlweDft trlweDft {param.k, param.N};

    genTrlweKey(trlweKey);

    // data gen
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }

    // enc
    Trglev trglev {param};
    encTrglevMultiSample(trglev, trlweKey, plainT, param);

    Integer y = genIntUniformDist(INT32_MIN, INT32_MAX);

    Trlwe recomp2 {param.k, param.N};
    multDecomposedTglevWithConst(recomp2, trglev, y, param);

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