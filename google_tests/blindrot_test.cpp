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

TEST(BlindRot, BlindRot) {
    YatfheParameters param {};
    param.n = 64;

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKey bsk {param};
    bootstrappingKeyGen(bsk, param, trgswKey, tlweKey);

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
    printArray(decIn.coeffs, "decIn");

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(INT32_MIN, INT32_MAX);
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
    trlweRotate(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    printArray(rotInP.coeffs, "expect");

    blindRotate(in2, bsk, sTlwe, param);

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

TEST(BlindRotLut, BlindRotLut) {
    YatfheParameters param {};
    param.n = 64;
    param.torusBase = 512;

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);
    BootstrappingKey bsk {param};
    bootstrappingKeyGen(bsk, param, trgswKey, tlweKey);
    TlweKeySwitchingKey ksk {param};
    tlweKeySwitchingKeyGen(ksk, trlweKey, tlweKey, param);

    // data gen
    Integer in = 3;
    Torus mu = modSwitchToTorus32(in, param.torusBase);
    cout << "in: " << in << endl;
    Tlwe tlwe {param.n};
    symEncTlweSample(tlwe, mu, tlweKey);
    auto decPre = symDecTlweSampleToInt(tlwe, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweFromTorus32(sTlwe, tlwe);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
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
    extractTlweFromTrlwe(tmp, in2, 0);
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