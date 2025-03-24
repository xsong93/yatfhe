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

TEST(BlindRot, BlindRot) {
    YatfheParameters param {};
    param.n = 64;
    yatfheInit(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKey bsk {param};
    bootstrappingKeyGen(bsk, trgswKey, tlweKey, param);

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
    cout << "rot:" << rot << endl;
    trlweRotate(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    blindRotate(in2, bsk, sTlwe, param);

    // trgsw enc X^rot
    Trgsw trgswXRot {param};
    Trlwe trlweInCopy {param.k, param.N};
    Trlwe trlwe {param.k, param.N};
    trgswEncrypt(trgswXRot, 1, trgswKey, 0, param);
    trgswRotate(trgswXRot, rot, param);
    symEncTrlweMultiSample(trlweInCopy, trlweKey, plainT.coeffs);
    trgswExternalProduct(trlwe, trgswXRot, trlweInCopy, param);
    IntPolynomial trlweDec {param.N};
    symDecTrlweToInt(trlweDec, trlwe, trlweKey, param.torusBase);
    printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    //verify
    for (auto i = 0; i < decP.N; i++) {
        ASSERT_EQ(rotInP.coeffs[i], decP.coeffs[i]);
    }
    printBanner("BlindRot");
}

TEST(BlindRot, BLIND_ROT_APPROX_CRT) {
    YatfheParameters param {};
    param.q = Q_CRT;
    param.n = 64;
    param.N = 32;
    yatfheInit(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKeyCRT bsKeyCRT{param};
    bootstrappingKeyGenApproxCRT(bsKeyCRT, trgswKey, tlweKey, param);

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
    trlweRotate(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    trlweMCRTDecomp(rotCRT, in2, param);
    blindRotateApproxCRT(rotCRT, bsKeyCRT, sTlwe, param);
    trlweMCRTToCRT(rotCRT, param);
    trlweCRTRecomp(resMCRT, rotCRT, param);

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

TEST(BlindRot, BLIND_ROT_APPROX_CRT_NTT) {
    YatfheParameters param {};
    param.q = Q_CRT;
//    param.n = 64;
//    param.N = 32;
    yatfheInit(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKeyCRT bsKeyCRT{param};
    bootstrappingKeyGenApproxCRT(bsKeyCRT, trgswKey, tlweKey, param);

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
    trlweRotate(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    trlweMCRTDecomp(rotCRT, in2, param);
    blindRotateApproxCRTNtt(rotCRT, bsKeyCRT, sTlwe, param);
    trlweMCRTToCRT(rotCRT, param);
    trlweCRTRecomp(resMCRT, rotCRT, param);

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

TEST(BlindRot, BlindRotLut) {
    YatfheParameters param {};
//    param.torusBase = 64;
    yatfheInit(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKey bsk {param};
    bootstrappingKeyGen(bsk, trgswKey, tlweKey, param);
    TlweKeySwitchingKey ksk {param};
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    tlweKeySwitchingKeyGen(ksk, trlweKey, tlweKsKey, param);

    // data gen
    Integer in = 3;
    Torus mu = modSwitchToTorus32(in, param.torusBase);
    cout << "in: " << in << endl;
    Tlwe tlwe {param.n};
    symEncTlweSample(tlwe, mu, tlweKey);
    auto decPre = symDecTlweSampleToInt(tlwe, tlweKey, param.torusBase);
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
    torusPolynomialRotate(rotIn, rot, in2.b);
    torusPolyToIntPoly(rotInP, rotIn, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    blindRotate(in2, bsk, sTlwe, param);
    printTrlweAB(in2, "aft brot");

    // trlwe dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    printArray(decP.coeffs, "real");

    // ks
    Tlwe tmp {ksk.nCurrKey};
    Tlwe tlweKs {ksk.nCurrKey};
    extractTlweFromTrlwe(tmp, in2, param.driftPhase);
    tlweKeySwitch(tlweKs, ksk, tmp, param);

    // tlwe dec
    auto out = symDecTlweSampleToInt(tlweKs, tlweKey, param.torusBase);
    cout << "out: " << out << endl;

    //verify
    for (auto i = 0; i < decP.N; i++) {
        ASSERT_NEAR(rotInP.coeffs[i], decP.coeffs[i], 1);
    }
    ASSERT_NEAR(in, out, 1);
    printBanner("BlindRotLut");
}