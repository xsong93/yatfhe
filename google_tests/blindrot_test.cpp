//
// Created by Xintong Song on 2024/4/29.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/blind_rotate.h"

TEST(BLIND_ROT, BLIND_ROT) {
    YatfheParameters param {};
    param.n = 3 * 20;
    param.group = 3;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param};
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
    // printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    // printArray(plain.coeffs, "plain");
    // printArray(decIn.coeffs, "decIn");

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
    // printMsg(rot, "rot");
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

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
    // printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    // printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT.BLIND_ROT");
}

TEST(BLIND_ROT, BLIND_ROT_NTT) {
    YatfheParameters param{};
    // param.n = 64;
    // param.group = 3;
    // param.n = 585;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param};
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
    COUNT_TIME("genBootstrappingKey", genBootstrappingKey(bsk, trgswKey, tlweKey, param);)

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
    // printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    // printArray(plain.coeffs, "plain");
    // printArray(decIn.coeffs, "decIn");

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
        }else if (tlweKey.s[i] == -1) {
            rot -= sTlwe.a[i];
        }
    }
    // printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

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
    // printArray(trlweDec.coeffs, "encXRot");

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    // printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT.BLIND_ROT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_MP21_NTT) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyMP bskMP{param, param.lApprox};
    COUNT_TIME("genBootstrappingKeyMP", genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);)

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    // pre dec
    IntPolynomial decIn{param.N};
    symDecTrlweToIntNtt(decIn, in2Dft, trlweKey, param.torusBase);
    // printArray(decIn.coeffs, "decIn");

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
        } else if (tlweKey.s[i] == -1) {
            rot -= sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

    COUNT_TIME("blindRotateMP21Ntt",
        blindRotateMP21Ntt(in2, bskMP.bskDft, sTlwe, param);)

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    // printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT.BLIND_ROT_MP21_NTT");
}

#ifdef TORUS32
static void benchBlindRotateMP21Ntt(int N, int radixBits, int l, int lApprox) {
    YatfheParameters param{};
    param.N = N;
    param.setRadixBits(radixBits);
    param.l = l;
    param.lApprox = lApprox;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);

    // data gen
    Trlwe in2{param.k, param.N};
    TrlweDft in2Dft{param.k, param.N};
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, v.coeffs);

    // rots gen
    ScaledTlwe sTlwe {param.N * 2, param.n};
    for (auto i = 0 ; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }

    // expected result
    int rot = 0;
    IntPolynomial rotInP{param.N};
    TrlweDft rotInDft{param.k, param.N};
    for (auto i = 0 ; i < param.n; i++) {
        if (tlweKey.s[i] == 1) {
            rot += sTlwe.a[i];
        } else if (tlweKey.s[i] == -1) {
            rot -= sTlwe.a[i];
        }
    }
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);

    std::string label = "blindRotateMP21Ntt N=" + std::to_string(param.N)
        + " radixBits=" + std::to_string(param.radixBits)
        + " l=" + std::to_string(param.l)
        + " lApprox=" + std::to_string(param.lApprox);
    // Wall-clock time (COUNT_TIME) plus total CPU work summed across all
    // threads via CLOCK_PROCESS_CPUTIME_ID, so it is independent of the
    // thread speed-up. cpu/wall approximates the effective parallelism.
    auto wallStart = steady_clock::now();
    timespec cpuStart{}, cpuEnd{};
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuStart);
    blindRotateMP21Ntt(in2, bskMP.bskDft, sTlwe, param);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuEnd);
    long wallUs = duration_cast<microseconds>(steady_clock::now() - wallStart).count();
    long cpuUs = (cpuEnd.tv_sec - cpuStart.tv_sec) * 1000000L
               + (cpuEnd.tv_nsec - cpuStart.tv_nsec) / 1000L;
    std::cout << "elapsed time (" << label << ") in us: " << wallUs << std::endl;
    std::cout << "total cpu work (" << label << ") in us: " << cpuUs
              << " (x" << (wallUs > 0 ? (double) cpuUs / wallUs : 0.0)
              << " parallelism)" << std::endl;

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);

    // soft correctness check: report mismatches but never fail (timing test)
    int mismatch = 0;
    for (auto i = 0; i < param.N; i++) {
        if (rotInP.coeffs[i] != decP.coeffs[i]) mismatch++;
    }
    if (mismatch == 0) {
        std::cout << "[correctness] " << label << ": OK" << std::endl;
    } else {
        std::cout << "[correctness] " << label << ": MISMATCH in "
                  << mismatch << "/" << param.N << " coeffs (not failing)"
                  << std::endl;
    }
}

// One test per (N, radixBits, l, lApprox) set so each runs a single
// initYatfhe (matching the rest of the suite) and can be filtered/timed
// independently.
TEST(BLIND_ROT, BLIND_ROT_MP21_NTT_TIMING_1) {
    benchBlindRotateMP21Ntt(2048, 7, 4, 3);
    printBanner("BLIND_ROT.BLIND_ROT_MP21_NTT_TIMING_1");
}

TEST(BLIND_ROT, BLIND_ROT_MP21_NTT_TIMING_2) {
    benchBlindRotateMP21Ntt(4096, 6, 5, 4);
    printBanner("BLIND_ROT.BLIND_ROT_MP21_NTT_TIMING_2");
}

TEST(BLIND_ROT, BLIND_ROT_MP21_NTT_TIMING_3) {
    benchBlindRotateMP21Ntt(8192, 4, 8, 6);
    printBanner("BLIND_ROT.BLIND_ROT_MP21_NTT_TIMING_3");
}

TEST(BLIND_ROT, BLIND_ROT_MP21_NTT_TIMING_4) {
    benchBlindRotateMP21Ntt(16384, 2, 16, 14);
    printBanner("BLIND_ROT.BLIND_ROT_MP21_NTT_TIMING_4");
}
#endif

TEST(BLIND_ROT, BLIND_ROT_PRE_ROT_NTT) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param};
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
    // printArray(decIn.coeffs, "decIn");

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
        } else if (tlweKey.s[i] == -1) {
            rot -= sTlwe.a[i];
        }
    }
    printMsg(rot, "rot");
    rotateTrlweNtt(rotInDft, in2Dft, rot);
    symDecTrlweToIntNtt(rotInP, rotInDft, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

    vector trgswMPDft(param.n, TrgswMPDft{param});
    COUNT_TIME("blindRotateWithPreRotNttMT",
               blindRotateWithPreRotNttMT(in2, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);)

    // dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    // printArray(decP.coeffs, "real");

    //verify
    ASSERT_EQ(rotInP.coeffs, decP.coeffs);
    printBanner("BLIND_ROT.BLIND_ROT_PRE_ROT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_APPROX_CRT) {
    YatfheParameters param {};
    param.q = Q_CRT;
    param.n = 64;
    param.N = 32;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n};
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
    // printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    // printArray(plain.coeffs, "plain");
    // printArray(decIn.coeffs, "decIn");

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
    // cout << "rot:" << rot << endl;
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

    decompTrlweMcrt(rotCRT, in2, param);
    blindRotateApproxCRT(rotCRT, bsKeyCRT.bsk8, sTlwe, param);
    trlweMcrtToCrt(rotCRT, param);
    recompTrlweCrt(resMCRT, rotCRT, param);

    // dec
    IntPolynomial decMP {param.N};
    symDecTrlweToInt(decMP, resMCRT, trlweKey, param.torusBase);
    // printArray(decMP.coeffs, "realMCRT");

    //verify
    for (auto i = 0; i < decMP.N; i++) {
        ASSERT_EQ(rotInP.coeffs[i], decMP.coeffs[i]);
    }
    printBanner("BLIND_ROT.BLIND_ROT_APPROX_CRT");
}

TEST(BLIND_ROT, BLIND_ROT_APPROX_CRT_NTT) {
    YatfheParameters param {};
    param.q = Q_CRT;
//    param.n = 64;
//    param.N = 32;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n};
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
    // printTrlweAB(in2, "input");

    // pre dec
    IntPolynomial decIn {param.N};
    symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
    // printArray(plain.coeffs, "plain");
    // printArray(decIn.coeffs, "decIn");

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
    // cout << "rot:" << rot << endl;
    rotateTrlwe(rotIn, in2, rot);
    symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
    // printArray(rotInP.coeffs, "expect");

    decompTrlweMcrt(rotCRT, in2, param);
    blindRotateApproxCRTNtt(rotCRT, bsKeyCRT.bskCRT, sTlwe, param);
    trlweMcrtToCrt(rotCRT, param);
    recompTrlweCrt(resMCRT, rotCRT, param);

    // dec
    IntPolynomial decMP {param.N};
    symDecTrlweToInt(decMP, resMCRT, trlweKey, param.torusBase);
    // printArray(decMP.coeffs, "realMCRT");

    //verify
    for (auto i = 0; i < decMP.N; i++) {
        ASSERT_EQ(rotInP.coeffs[i], decMP.coeffs[i]);
    }
    printBanner("BLIND_ROT.BLIND_ROT_APPROX_CRT_NTT");
}

TEST(BLIND_ROT, BLIND_ROT_LUT) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyMP bsk {param};
    genBootstrappingKeyMP(bsk, trgswKey, tlweKey, param);
    TlweKeySwitchingKey ksk {param};
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksk, trlweKey, tlweKsKey, param);

    // data gen
    Integer in = 3;
    Torus mu = modSwitchToTorusGeneral(in, param.torusBase, LWE_Q);
    // cout << "in: " << in << endl;
    Tlwe tlwe {param.n};
    symEncTlwe(tlwe, mu, tlweKey);
    auto decPre = symDecTlweToInt(tlwe, tlweKey, param.torusBase);
    // cout << "decPre: " << decPre << endl;
    ScaledTlwe sTlwe {param.N * 2, param.n};
    // printTlweAB(tlwe, "tlwe");
    rescaleTlweToNewMod(sTlwe, tlwe);
    // printTlweAB(sTlwe, "sTlwe");
    // printArray(tlweKey.s, "s");

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    // printArray(v.coeffs, "v");
    Trlwe in2 {param.k, param.N};
    genNoiselessTrlweSample(in2, v, sTlwe);
    // printTrlweAB(in2, "input");

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
    // printArray(rotInP.coeffs, "expect");

    COUNT_TIME("blindRotate", blindRotateMP21Ntt(in2, bsk.bskDft, sTlwe, param);)
    // printTrlweAB(in2, "aft brot");

    // trlwe dec
    IntPolynomial decP {param.N};
    symDecTrlweToInt(decP, in2, trlweKey, param.torusBase);
    // printArray(decP.coeffs, "real");

    // ks
    Tlwe tmp {ksk.nCurrKey};
    Tlwe tlweKs {ksk.nCurrKey};
    extractTlweFromTrlwe(tmp, in2, param.driftPhase);
    switchKeyForTlwe(tlweKs, ksk, tmp, param);

    // tlwe dec
    auto out = symDecTlweToInt(tlweKs, tlweKey, param.torusBase);
    // cout << "out: " << out << endl;

    //verify
    for (auto i = 0; i < decP.N; i++) {
        ASSERT_NEAR(rotInP.coeffs[i], decP.coeffs[i], 1);
    }
    ASSERT_NEAR(in, out, 1);
    printBanner("BLIND_ROT.BLIND_ROT_LUT");
}

TEST(BLIND_ROT, V2IR_MP21_NTT) {
    YatfheParameters param{};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyMP bskMP{param, param.lApprox};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);

    int pt;
    for (pt = -param.torusBase/2; pt < param.torusBase/2; pt++) {
        // data gen
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);

        // pre dec
        auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
        cout << "decPre: " << decPre << endl;

        // v2i
        Trlwe output{param.k, param.N};
        COUNT_TIME("transferValueToIndexRange", transferValueToIndexRange(output, input, 1, bskMP, param);)

        // dec
        IntPolynomial decP{param.N};
        symDecTrlweToInt(decP, output, trlweKey, param.torusBase);
        // printArray(decP.coeffs, "real");
        int count = 0;
        cout << "non zero: ";
        for (auto i = 0; i < param.N; i++) {
            if (decP.coeffs[i] != 0) {
                cout << i << ":" << decP.coeffs[i] << ", ";
                count++;
            }
        }
        cout << endl << "non zero count: " << count << endl;

        //verify
        auto scale = 2 * param.N / param.torusBase;
        // cout << "search index:" << param.N - 1 - scale * pt << endl;
        // ASSERT_EQ(decP.coeffs[param.N - 1 - scale * pt], -1);
        int sign = pt < 0 ? -1 : 1;
        int searchIdx =  pt < 0 ? -scale * pt - 1 : param.N - 1 - scale * pt;
        cout << "searchIdx: " << searchIdx << endl;
        ASSERT_EQ(decP.coeffs[searchIdx], -1 * sign);
        ASSERT_EQ(count, scale);
    }
    printBanner("BLIND_ROT.V2IR_MP21_NTT");
}