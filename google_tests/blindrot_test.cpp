//
// Created by Xintong Song on 2024/4/29.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/ntt_hexl.h"

TEST(BLIND_ROT, BLIND_ROT) {
    YatfheParameters param {};
    param.n = 3 * 20;
    param.group = 3;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKey bsk {param};
    genBootstrappingKey(bsk, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2 {param.k, param.N};
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }
    symEncTrlweMultiSample(in2, trlweKey, plainT.coeffs);
    printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    printArray(plain.coeffs, "plain");
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    }

    // test data gen
    int rot = 0;
    IntPolynomial rotInP {param.N};
    Trlwe rotIn {param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    COUNT_TIME("blindRotate", blindRotate(in2, bsk.bsk, sTlwe, param);)

    // trgsw enc X^rot
    Trgsw trgswXRot {param};
    Trlwe trlweInCopy {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    encryptTrgsw(trgswXRot, 1, trgswKey, 0, param);
    rotateTrgsw(trgswXRot, rot, param);
    symEncTrlweMultiSample(trlweInCopy, trlweKey, plainT.coeffs);
    externalProductTrgsw(trlwe, trgswXRot, trlweInCopy, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT");
}

TEST(BLIND_ROT, BLIND_ROT_NTT) {
    YatfheParameters param{};
    // param.n = 64;
    // param.group = 3;
    // param.n = 585;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
#ifdef DEBUG_MODE
    for (int i = 0; i < tlweKey.n; i = i + 2) {
        tlweKey.s[i] = 0;
        tlweKey.s[i + 1] = 1;
    }
#endif
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKey bsk {param};
    genBootstrappingKey(bsk, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    IntPolynomial plain{param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, plainT.coeffs);
    printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    printArray(plain.coeffs, "plain");
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // test data gen
    int rot = 0;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

#ifdef DEBUG_MODE
    IntPolynomial ip1{param.N}, ip2{param.N};
    Trlwe temp{param.k, param.N}, temp2{param.k, param.N};
    TrgswDft tmp1{param}, tmp2{param}, tmp3{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    for (auto i = 0; i < param.n; i = i + 2) {
        auto a1 = sTlwe.a[i];
        auto a2 = sTlwe.a[i + 1];
        auto bsk1 = bsk.bskDft[j];
        auto& bsk2 = bsk.bskDft[j + 1];
        auto bsk3 = bsk.bskDft[j + 2];
        auto& bsk4 = bsk.bskDft[j + 3];

        rotateTrgswNtt(bsk1, a1, param);
        rotateTrgswNtt(bsk3, a1, param);
        addTrgswNtt(tmp1, bsk1, bsk2);
        addTrgswNtt(tmp2, bsk3, bsk4);

        rotateTrgswNtt(tmp1, a2, param);
        addTrgswNtt(tmp3, tmp1, tmp2);

        temp = Trlwe{param.k, param.N};
        externalProductTrgswNtt(temp, tmp3, in2, param);
        rotateTrlwe(temp2, in2, a2);
        symDecTrlweToInt(ip1, temp, trlweKey, param.torusBase);
        symDecTrlweToInt(ip2, temp2, trlweKey, param.torusBase);
        try {
            if (ip1.coeffs != ip2.coeffs) throw std::runtime_error("");
        } catch (...) {
            printMsg(a2, "a2");
            ASSERT_EQ(ip1.coeffs, ip2.coeffs);
        }
        in2 = std::move(temp);
        j += batchSize;
    }
#else
    COUNT_TIME("blindRotateNtt", blindRotateNtt(in2, bsk.bskDft, sTlwe, param);)
#endif

    // trgsw enc X^rot
    Trgsw trgswXRot{param};
    TrgswDft trgswXRotDft{param};
    Trlwe trlweInCopy{param.k, param.N};
    TrlweDft trlweInCopyDft{param.k, param.N};
    Trlwe trlwe{param.k, param.N};
    encryptTrgswNtt(trgswXRot, trgswXRotDft, 1, trgswKey, 0, param);
    rotateTrgswNtt(trgswXRotDft, rot, param);
    symEncTrlweMultiSampleNtt(trlweInCopy, trlweInCopyDft, trlweKey, plainT.coeffs);
    externalProductTrgswNtt(trlwe, trgswXRotDft, trlweInCopy, param.lApprox, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_PRE_ROT_NTT) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    BootstrappingKeyMPPreRot bskPre{param};
    genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    sTlwe.b = genIntUniformDist(-2 * param.N, 2 * param.N);
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // test data gen
    int rot = -sTlwe.b;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    vector trgswMPDft(param.n, TrgswMPDft{param});
    COUNT_TIME("blindRotateWithPreRotNtt",
               blindRotateWithPreRotNtt(in2, bskPre.bskFirst, bskPre.bskDft, sTlwe, param.batchSize, param);)

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT_PRE_ROT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_INTERNAL_PAIRWISE_ASYM_OPT_NTT) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    BootstrappingKeyInternalAsymOpt bskAsymOpt{param};
    genBootstrappingKeyInternalAsymOpt(bskAsymOpt, trgswKey, tlweKey, v, param);

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    sTlwe.b = 0;
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // test data gen
    int rot = -sTlwe.b;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    TrgswMP accDummy{param};
    encryptTrgswMP(accDummy, 1, trgswKey, 0, param);
    vector bskDummy(param.n, TrgswMPDft(param));
    for (auto i = 0; i < param.n; i++) {
        encryptTrgswMPNtt(bskDummy[i], 1, trgswKey, 0, param);
    }
    TrlevDft s2Dft(param);
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

    COUNT_TIME("blindRotateInternalPairWiseAsymOptNtt",
               blindRotateInternalPairWiseAsymOptNtt(in2, bskAsymOpt.bsk, bskAsymOpt.bskLast, bskAsymOpt.bskDft,
                                                     sTlwe, s2Dft, param.batchSize, param);)

    // trgsw enc X^rot
    Trgsw trgswXRot{param};
    TrgswDft trgswXRotDft{param};
    Trlwe trlweInCopy{param.k, param.N};
    TrlweDft trlweInCopyDft{param.k, param.N};
    Trlwe trlwe{param.k, param.N};
    encryptTrgswNtt(trgswXRot, trgswXRotDft, 1, trgswKey, 0, param);
    rotateTrgswNtt(trgswXRotDft, rot, param);
    symEncTrlweMultiSampleNtt(trlweInCopy, trlweInCopyDft, trlweKey, v.coeffs);
    externalProductTrgswNtt(trlwe, trgswXRotDft, trlweInCopy, param.lApprox, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT_INTERNAL_PAIRWISE_ASYM_OPT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_INTERNAL_PAIRWISE_ASYM_NTT) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
//    BootstrappingKeyInternal bsk{param};
//    genBootstrappingKeyInternal(bsk, trgswKey, tlweKey, param);
//    BootstrappingKey bskNor{param};
//    genBootstrappingKey(bskNor, trgswKey, tlweKey, param);
    BootstrappingKeyInternalAsym bskAsym{param};
    genBootstrappingKeyInternalAsym(bskAsym, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
//    IntPolynomial plain{param.N}; // Z/pZ
//    TorusPolynomial plainT {param.N};
//    for (auto i = 0; i < plain.N; i++) {
//        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
//        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
//    }
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    sTlwe.b = 0;
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // test data gen
    int rot = -sTlwe.b;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    TrgswMP accDummy{param};
    encryptTrgswMP(accDummy, 1, trgswKey, 0, param);
    vector bskDummy(param.n, TrgswMPDft(param));
    for (auto i = 0; i < param.n; i++) {
        encryptTrgswMPNtt(bskDummy[i], 1, trgswKey, 0, param);
    }
    TrlevDft s2Dft(param);
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

//    COUNT_TIME("blindRotateNtt", blindRotateNtt(in2, bskNor.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateInternalNtt", blindRotateMPInternalNtt(accDummy, bskDummy, sTlwe, param);)
//    COUNT_TIME("blindRotateInternalPireWiseNtt", blindRotateInternalPairWiseNtt(in2, bsk.bsk, bsk.bskDft, sTlwe, param);)
    COUNT_TIME("blindRotateInternalPairWiseAsymNtt", blindRotateInternalPairWiseAsymNtt(in2, bskAsym.bsk, bskAsym.bskDft, sTlwe, s2Dft, param.batchSize, param);)

    // trgsw enc X^rot
    Trgsw trgswXRot{param};
    TrgswDft trgswXRotDft{param};
    Trlwe trlweInCopy{param.k, param.N};
    TrlweDft trlweInCopyDft{param.k, param.N};
    Trlwe trlwe{param.k, param.N};
    encryptTrgswNtt(trgswXRot, trgswXRotDft, 1, trgswKey, 0, param);
    rotateTrgswNtt(trgswXRotDft, rot, param);
    symEncTrlweMultiSampleNtt(trlweInCopy, trlweInCopyDft, trlweKey, v.coeffs);
    externalProductTrgswNtt(trlwe, trgswXRotDft, trlweInCopy, param.lApprox, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT_INTERNAL_PAIRWISE_ASYM_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_INTERNAL_PAIRWISE_NTT) {
    YatfheParameters param{};
//    param.n = 4;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyInternal bsk{param};
    genBootstrappingKeyInternal(bsk, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    sTlwe.b = 0;
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // test data gen
    int rot = -sTlwe.b;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

//    COUNT_TIME("blindRotateInternalNtt", blindRotateMPInternalNtt(accDummy, bskDummy, sTlwe, param);)
    COUNT_TIME("blindRotateInternalPireWiseNtt",
               blindRotateInternalPairWiseNtt(in2, bsk.bsk, bsk.bskDft, sTlwe, param.batchSize, param);)

    // trgsw enc X^rot
    Trgsw trgswXRot{param};
    TrgswDft trgswXRotDft{param};
    Trlwe trlweInCopy{param.k, param.N};
    TrlweDft trlweInCopyDft{param.k, param.N};
    Trlwe trlwe{param.k, param.N};
    encryptTrgswNtt(trgswXRot, trgswXRotDft, 1, trgswKey, 0, param);
    rotateTrgswNtt(trgswXRotDft, rot, param);
    symEncTrlweMultiSampleNtt(trlweInCopy, trlweInCopyDft, trlweKey, v.coeffs);
    externalProductTrgswNtt(trlwe, trgswXRotDft, trlweInCopy, param.lApprox, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT_INTERNAL_PAIRWISE_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_APPROX_CRT) {
    YatfheParameters param {};
    param.q = Q_CRT;
    param.n = 64;
    param.N = 32;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyCRT bsKeyCRT{param};
    genBootstrappingKeyApproxCrt(bsKeyCRT, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2 {param.k, param.N};
    Trlwe resMCRT {param.k, param.N};
    std::vector<Trlwe8> rotCRT(param.d, Trlwe8{param.k, param.N});
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }
    symEncTrlweMultiSample(in2, trlweKey, plainT.coeffs);
    printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    printArray(plain.coeffs, "plain");
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    }

    // test data gen
    int rot = 0;
    IntPolynomial rotInP {param.N};
    Trlwe rotIn {param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    cout << "rot:" << rot << endl;
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    decompTrlweMcrt(rotCRT, in2, param);
    blindRotateApproxCRT(rotCRT, bsKeyCRT.bsk8, sTlwe, param);
    trlweMcrtToCrt(rotCRT, param);
    recompTrlweCrt(resMCRT, rotCRT, param);

    // dec
    IntPolynomial decMP {param.N};
    symDecTrlweToInt(decMP, resMCRT, trlweKey, param.torusBase);
    printArray(decMP.coeffs, "realMCRT");

    //verify
    for (auto i = 0; i < decMP.N; i++) {
        ASSERT_EQ(rotInP.coeffs[i], decMP.coeffs[i]);
    }
    printBanner("BLIND_ROT_APPROX_CRT");
}

TEST(BLIND_ROT, BLIND_ROT_APPROX_CRT_NTT) {
    YatfheParameters param {};
    param.q = Q_CRT;
//    param.n = 64;
//    param.N = 32;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyCRT bsKeyCRT{param};
    genBootstrappingKeyApproxCrt(bsKeyCRT, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2 {param.k, param.N};
    Trlwe resMCRT {param.k, param.N};
    std::vector<Trlwe8> rotCRT(param.d, Trlwe8{param.k, param.N});
    IntPolynomial plain {param.N}; // Z/pZ
    TorusPolynomial plainT {param.N};
    for (auto i = 0; i < plain.N; i++) {
        plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
        plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
    }
    symEncTrlweMultiSample(in2, trlweKey, plainT.coeffs);
    printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    printArray(plain.coeffs, "plain");
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    }

    // test data gen
    int rot = 0;
    IntPolynomial rotInP {param.N};
    Trlwe rotIn {param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    cout << "rot:" << rot << endl;
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    decompTrlweMcrt(rotCRT, in2, param);
    blindRotateApproxCRTNtt(rotCRT, bsKeyCRT.bskCRT, sTlwe, param);
    trlweMcrtToCrt(rotCRT, param);
    recompTrlweCrt(resMCRT, rotCRT, param);

    // dec
    IntPolynomial decMP {param.N};
    symDecTrlweToInt(decMP, resMCRT, trlweKey, param.torusBase);
    printArray(decMP.coeffs, "realMCRT");

    //verify
    for (auto i = 0; i < decMP.N; i++) {
        ASSERT_EQ(rotInP.coeffs[i], decMP.coeffs[i]);
    }
    printBanner("BLIND_ROT_APPROX_CRT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_LUT) {
    YatfheParameters param {};
    param.n = 64;
    param.group = 2;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKey bsk {param};
    genBootstrappingKey(bsk, trgswKey, tlweKey, param);
    TlweKeySwitchingKey ksk {param};
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksk, trlweKey, tlweKsKey, param);

    // data gen
    Integer in = 3;
    Torus mu = modSwitchToTorus32(in, param.torusBase);
    cout << "in: " << in << endl;
    Tlwe tlwe {param.n};
    symEncTlwe(tlwe, mu, tlweKey);
    auto decPre = symDecTlweToInt(tlwe, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;
    ScaledTlwe sTlwe {param.N * 2, param.n};
    printTlweAB(tlwe, "tlwe");
    rescaleTlweFromTorus32(sTlwe, tlwe);
    printTlweAB(sTlwe, "sTlwe");
    printArray(tlweKey.s, "s");

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    printArray(v.coeffs, "v");
    Trlwe in2 {param.k, param.N};
    genNoiselessTrlweSample(in2, v, sTlwe);
    printTrlweAB(in2, "input");

    // test data gen
    int rot = 0;
    IntPolynomial rotInP {param.N};
    TorusPolynomial rotIn {param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        }
    }
    rotateTorusPolynomial(rotIn, rot, in2.b);
    torusPolyToIntPoly(rotInP, rotIn, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    COUNT_TIME("blindRotate", blindRotate(in2, bsk.bsk, sTlwe, param);)
    printTrlweAB(in2, "aft brot");

    // trlwe dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    // ks
    Tlwe tmp {ksk.nCurrKey};
    Tlwe tlweKs {ksk.nCurrKey};
    extractTlweFromTrlwe(tmp, in2, param.driftPhase);
    switchKeyForTlwe(tlweKs, ksk, tmp, param);

    // tlwe dec
    auto out = symDecTlweToInt(tlweKs, tlweKey, param.torusBase);
    cout << "out: " << out << endl;

    //verify
    for (auto i = 0; i < decP.N; i++) {
        ASSERT_NEAR(rotInP.coeffs[i], decP.coeffs[i], 1);
    }
    ASSERT_NEAR(in, out, 1);
    printBanner("BLIND_ROT_LUT");
}