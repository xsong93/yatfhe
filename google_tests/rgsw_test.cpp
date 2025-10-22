//
// Created by Xintong Song on 2024/3/26.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlgsw.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/ntt_hexl.h"

TEST(RgswTest, RgswEncDecTest) {
    YatfheParameters param {};
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain = 7;
    encryptTrgswNtt(trgsw, trgswDft, plain, trgswKey, 0, param);
//    trgswEncrypt(trgsw, param, trgswKey, plain);

//
//    // test identity for trgsw and trgswDft value
//    for (auto i = 0; i < trgsw.l; i++) {
//        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
//            TorusPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
//            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
//            ASSERT_EQ(ip.coeffs, trgsw.trlweSamples[i][j].b.coeffs);
//        }
//    }

    // trgsw dec
    Integer dec = decryptTrgswNtt(trgswDft, param, trgswKey);
//    Integer dec = trgswDecrypt(trgsw, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RgswEncDecTest");
}

TEST(RgswTest, RGSW_ADD) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw1{param}, trgsw2{param}, res{param};
    TrgswDft trgswDft1{param}, trgswDft2{param}, resDft{param};
    Integer dec, plainAdd;
    Integer plain1, plain2;
    plain1 = 1;
    plain2 = 0;
    encryptTrgswNtt(trgsw1, trgswDft1, plain1, trgswKey, 0, param);
    encryptTrgswNtt(trgsw2, trgswDft2, plain2, trgswKey, 0, param);
    addTrgsw(res, trgsw1, trgsw2);
    addTrgswNtt(resDft, trgswDft1, trgswDft2);

    // trgsw dec
    dec = decryptTrgsw(res, param, trgswKey);
    plainAdd = intModP(plain1 + plain2, param.torusBase);
    cout << "plain: " << plainAdd << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plainAdd, dec);

    dec = decryptTrgswNtt(resDft, param, trgswKey);
    ASSERT_EQ(plainAdd, dec);

    plain1 = 3;
    plain2 = -1;
    encryptTrgswNtt(trgsw1, trgswDft1, plain1, trgswKey, 0, param);
    encryptTrgswNtt(trgsw2, trgswDft2, plain2, trgswKey, 0, param);
    addTrgsw(res, trgsw1, trgsw2);
    addTrgswNtt(resDft, trgswDft1, trgswDft2);

    // trgsw dec
    dec = decryptTrgsw(res, param, trgswKey);
    plainAdd = intModP(plain1 + plain2, param.torusBase);
    cout << "plain: " << plainAdd << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plainAdd, dec);

    dec = decryptTrgswNtt(resDft, param, trgswKey);
    ASSERT_EQ(plainAdd, dec);
    printBanner("RGSW_ADD");
}

TEST(RgswTest, RGSW_SUB) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw1{param}, trgsw2{param}, res{param};
    TrgswDft trgswDft1{param}, trgswDft2{param}, resDft{param};
    Integer dec, plainSub;
    Integer plain1, plain2;
    plain1 = 1;
    plain2 = 0;
    encryptTrgswNtt(trgsw1, trgswDft1, plain1, trgswKey, 0, param);
    encryptTrgswNtt(trgsw2, trgswDft2, plain2, trgswKey, 0, param);
    subTrgsw(res, trgsw1, trgsw2);
    subTrgswNtt(resDft, trgswDft1, trgswDft2);

    // trgsw dec
    dec = decryptTrgsw(res, param, trgswKey);
    plainSub = intModP(plain1 - plain2, param.torusBase);
    cout << "plain: " << plainSub << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plainSub, dec);

    dec = decryptTrgswNtt(resDft, param, trgswKey);
    ASSERT_EQ(plainSub, dec);

    plain1 = 1;
    plain2 = -1;
    encryptTrgswNtt(trgsw1, trgswDft1, plain1, trgswKey, 0, param);
    encryptTrgswNtt(trgsw2, trgswDft2, plain2, trgswKey, 0, param);
    subTrgsw(res, trgsw1, trgsw2);
    subTrgswNtt(resDft, trgswDft1, trgswDft2);

    // trgsw dec
    dec = decryptTrgsw(res, param, trgswKey);
    plainSub = intModP(plain1 - plain2, param.torusBase);
    cout << "plain: " << plainSub << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plainSub, dec);

    dec = decryptTrgswNtt(resDft, param, trgswKey);
    ASSERT_EQ(plainSub, dec);
    printBanner("RGSW_SUB");
}

TEST(RgswTest, RGSW_MCRT_DECOMPOSITION) {
    YatfheParameters param {};
    param.q = Q_CRT;
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain = 7;
    encryptTrgsw(trgsw, plain, trgswKey, 0, param);

    // GD
    std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
    Trgsw trgswRe {param};
    decompTrgswMcrt(trgswD, trgsw, param);
    trgswMcrtToCrt(trgswD, param);
    recompTrgswCrt(trgswRe, trgswD, param);
    printTrgsw(trgsw, "trgsw");
    printTrgsw(trgswRe, "trgswRe");

    // trgsw dec
    Integer dec = decryptTrgsw(trgswRe, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RGSW_MCRT_DECOMPOSITION");
}

TEST(RgswTest, RgswMultTestNaive) {
    YatfheParameters param {};
//    param.N = 1024;
//    param.radixBits = 4;
//    param.l = 3;
//    param.k = 1;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        encryptTrgsw(trgsw, mu1, trgswKey, 0, param);
        printf( "trgsw dec: %d.\n", decryptTrgsw(trgsw, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = intModP(i, param.torusBase);
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu1, param.torusBase);
        }
        Trlwe out {param.k, param.N};
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult
        COUNT_TIME("trgswExternalProduct", externalProductTrgsw(out, trgsw, in2, param);)
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNaive");
}

TEST(RgswTest, RgswMultTestNTT) {
    YatfheParameters param {};
//    param.k = 8;
//    param.N = 128;
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
        Integer mu1 = genIntUniformDist(0, 3);
        encryptTrgswNtt(trgsw, trgswDft, mu1, trgswKey, 0, param);
        printf( "trgsw dec: %d.\n", decryptTrgswNtt(trgswDft, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft {param.k, param.N};
//        Integer mu2p = genIntUniformDist(0, 3);
//        Torus mu2 = modSwitchToTorus32(mu2p, param.torusBase);
        vector<Integer> mu2ps(param.N);
        vector<Torus> mu2s(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2ps[i] = genIntUniformDist(0, 3);
            mu2s[i] = modSwitchToTorus32(mu2ps[i], param.torusBase);
        }
        Trlwe out {param.k, param.N};
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
//        symEncTrlweSingleSampleNtt(in2, in2Dft, trlweKey, mu2);
        symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, mu2s);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        COUNT_TIME("trgswExternalProductNtt", externalProductTrgswNtt(out, trgswDft, in2, param.lApprox, param);)
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");

        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNTT");
}

TEST(RgswTest, RGSWMP_MULT_NAIVE) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
//    param.l = 2;
//    param.k = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        TrgswMP trgswMP {param};
        Integer mu1 = 1;
        encryptTrgswMP(trgswMP, mu1, trgswKey, 0, param);

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu1, param.torusBase);
        }
        Trlwe out {param.k, param.N};
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
//        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        COUNT_TIME("trgswMPExternalProduct", externalProductTrgswMP(out, trgswMP, in2, param.lApprox, param);)
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_MULT_NAIVE");
}

TEST(RgswTest, RGSWMP_MULT_NTT) {
    YatfheParameters param {};
//    param.k = 8;
//    param.N = 128;
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        TrgswMPDft trgswDft {param};
        Integer mu1 = genIntUniformDist(0, 3);
        encryptTrgswMPNtt(trgswDft, mu1, trgswKey, 0, param);

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft {param.k, param.N};
        vector<Integer> mu2ps(param.N);
        vector<Torus> mu2s(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2ps[i] = genIntUniformDist(0, 3);
            mu2s[i] = modSwitchToTorus32(mu2ps[i], param.torusBase);
        }
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, mu2s);
        printTrlweAB(in2, "trlwe");
        Trlwe out {param};

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        COUNT_TIME("trgswMPExternalProductNtt", externalProductTrgswMPNtt(out, trgswDft, in2, param.lApprox, param);)

        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");

        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_MULT_NTT");
}

TEST(RgswTest, RGSWMP_INTERMULT_NAIVE) {
    YatfheParameters param {};
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        TrgswMP trgswMP1 {param};
        Integer mu1 = 2;
        encryptTrgswMP(trgswMP1, mu1, trgswKey, 0, param);

        TrgswMP trgswMP2 {param};
        Integer mu2 = 1;
        encryptTrgswMP(trgswMP2, mu2, trgswKey, 0, param);

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu1, param.torusBase);
        }
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
        //        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        TrgswMP tmp{param};
        Trlwe out{param.k, param.N};
        TorusPolynomial t1{param.N};
        decryptTrgswMP(t1, trgswMP1, param, trgswKey, false);
        printArray(t1.coeffs, "trgswMP1");
        decryptTrgswMP(t1, trgswMP2, param, trgswKey, false);
        printArray(t1.coeffs, "trgswMP2");
        for (auto i = 0; i < 1; i++) {
            COUNT_TIME("trgswMPInternalProduct", internalProductTrgswMP(trgswMP1, trgswMP1, trgswMP2, param.l, param);)
        }
        COUNT_TIME("trgswMPExternalProduct", externalProductTrgswMP(out, trgswMP1, in2, param.lApprox, param);)


        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_INTERMULT_NAIVE");
}

TEST(RgswTest, RGSWMP_INTERMULT_NAIVE_NTT) {
    YatfheParameters param {};
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        TrgswMP trgswMP1{param};
        Integer mu1 = 2;
        encryptTrgswMP(trgswMP1, mu1, trgswKey, 0, param);

        TrgswMPDft trgswMP2{param};
        Integer mu2 = 1;
        encryptTrgswMPNtt(trgswMP2, mu2, trgswKey, 0, param);

        // trlwe enc
        Trlwe in2{param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = 1;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu1 * mu2, param.torusBase);
        }
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
        //        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP{param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        TrgswMP tmp{param};
        Trlwe out{param.k, param.N};
        TorusPolynomial t1{param.N};
        decryptTrgswMP(t1, trgswMP1, param, trgswKey, false);
        printArray(t1.coeffs, "trgswMP1");
        decryptTrgswMPNtt(t1, trgswMP2, param, trgswKey, false);
        printArray(t1.coeffs, "trgswMP2");
        int test = mu1;
        for (auto i = 0; i < 8; i++) {
            COUNT_TIME("internalProductTrgswMPNtt", internalProductTrgswMPNtt(trgswMP1, trgswMP1, trgswMP2, param.l, param);)
            test *= mu2;
        }
        cout << test << endl;
        COUNT_TIME("trgswMPExternalProduct", externalProductTrgswMP(out, trgswMP1, in2, param.lApprox, param);)


        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_INTERMULT_NAIVE_NTT");
}

TEST(RgswTest, RGSWMP_SCHEME_SWITCHING) {
    YatfheParameters param {};
    param.torusBits = 32;
    initYatfhe(param);
    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // enc s^2
    Trlev s2p(param, param.l);
    TrlevDft s2pDft(param, param.l);
    symEncTrlevWithKey(s2p, trlweKey, trlweKey.s, true, param);
    symEncTrlevWithKeyNtt(s2pDft, trlweKey, trlweKey.s, true, param);

    // trgsw enc
    TrgswMPDft in1Dft{param, param.lApprox};
    TrgswMP in1{param, param.lApprox, true};
    Integer mu1 = 1;
    encryptTrgswMP(in1, mu1, trgswKey, 0, param);
    encryptTrgswMPNtt(in1Dft, mu1, trgswKey, 0, param);

    IntPolynomial dec{param.N};
    decryptTrgswMPNtt(dec, in1Dft, param, trgswKey, false);
    printArray(dec.coeffs, "mu1");
    decryptTrgswMPNtt(dec, in1Dft, param, trgswKey, true);
    printArray(dec.coeffs, "-mu1*s");

    // init data
    TrgswMPDft tmpMp{param, param.lApprox};
    tmpMp.cPrime = in1Dft.cPrime;
    TrgswMPDft tmpMp2{param, param.lApprox};
    vector decompA(param.lApprox, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    for(auto l = 0; l < param.lApprox; l++) {
        NttHexl::applyNtt(tmpMp2.cPrime[l].b, in1.cPrime[l].b);
    }
    for (auto l0 = 0; l0 < param.lApprox; l0++) {
        auto& a = in1.cPrime[l0].a;
        for (auto k = 0; k < param.k; k++) {
            for (auto j = 0; j < param.N; j++) {
                DecomposedData d{param.l};
                gadgetDecompose(d, a[k].coeffs[j], param);
                for (auto l = 0; l < param.l; l++) {
                    decompA[l0][l][k].coeffs[j] = d.value[l] * d.sign;
                }
            }
        }
    }
    TrgswMP t1{param, param.lApprox};
    t1.cPrime = in1.cPrime;

    // scheme switching
    for (auto l = 0; l < param.lApprox; l++) {
        COUNT_TIME("switchTrlweToSecretEmbeddingNtt",
                   switchTrlweToSecretEmbeddingNtt(tmpMp.c[l], tmpMp.cPrime[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNttOpt",
                   switchTrlweToSecretEmbeddingNttOpt(tmpMp2.c[l], tmpMp2.cPrime[l], decompA[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbedding", switchTrlweToSecretEmbedding(t1.c[l], in1.cPrime[l], s2p, param);)
    }

    // result validation
    Trlwe trlwe{param.k, param.N};
    Trlwe res{param.k, param.N};
    TorusPolynomial mu{param.N};
    IntPolynomial plainMult{param.N};
    for (auto j = 0; j < param.N; j++) {
        plainMult.coeffs[j] = intModP(j * mu1, param.torusBase);
        mu.coeffs[j] = modSwitchToTorus32(j, param.torusBase);
    }
    symEncTrlweMultiSample(trlwe, trlweKey, mu.coeffs);

    // switchTrlweToSecretEmbedding
    externalProductTrgswMP(res, t1, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    ASSERT_EQ(dec.coeffs, plainMult.coeffs);

    // switchTrlweToSecretEmbeddingNtt
    externalProductTrgswMPNtt(res, tmpMp, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    ASSERT_EQ(dec.coeffs, plainMult.coeffs);

    // switchTrlweToSecretEmbeddingNttOpt
    externalProductTrgswMPNtt(res, tmpMp2, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    ASSERT_EQ(dec.coeffs, plainMult.coeffs);

    printBanner("RGSWMP_SCHEME_SWITCHING");
}

TEST(RgswTest, RGSWMP_INTERMULT_NTT) {
    YatfheParameters param {};
    param.torusBits = 32;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);
        TrlevDft s2Dft(param);
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

        // trgsw enc
        TrgswMPDft in1Dft{param};
        Integer mu1 = 3;
        encryptTrgswMPNtt(in1Dft, mu1, trgswKey, 0, param);
        IntPolynomial dec{param.N};
        decryptTrgswMPNtt(dec, in1Dft, param, trgswKey, false);
        printArray(dec.coeffs, "mu1");

        // trlwe enc
        Trlev in2 {param};
        Trlwe in3{param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        TorusPolynomial mu2T{param.N};
        for (auto i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i == 0 ? 1 : 0;
            multPlain.coeffs[i] = modMulQ(i, mu1, param.torusBase);
            mu2T.coeffs[i] = modSwitchToTorus32(i, param.torusBase);
        }
        encTrlevMultiSample(in2, trlweKey, mu2p, param);
        decTrlev(dec, in2, trlweKey, param);
        printArray(dec.coeffs, "in2");
        symEncTrlweMultiSample(in3, trlweKey, mu2T.coeffs);

        // trgsw mult
        TrgswMPDft tmp{param};
        COUNT_TIME("internalProductAsymTrgswMPNtt", internalProductAsymTrgswMPNtt(tmp, in1Dft, in2, s2Dft, param.lApprox, param))
        Trlwe res{param.k, param.N};
        externalProductTrgswMPNtt(res, tmp, in3, param.lApprox, param);

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        decryptTrgswMPNtt(dec, tmp, param, trgswKey, false);
        symDecTrlweToInt(decAftP, res, trlweKey, param.torusBase);
        printArray(dec.coeffs, "dec");
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        ASSERT_EQ(decAftP.coeffs, multPlain.coeffs);
    }
    printBanner("RGSWMP_INTERMULT_NTT");
}

// max depth 4
TEST(RgswTest, RGSW_MULT_NAIVE_CHAIN_NTT) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
    // param.l = 4;
    // param.lweStdDev = 0;
    // param.rlweStdDev = 0;
//    param.k = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        int loop = 4;
        std::vector<Trgsw> trgsws(loop, Trgsw{param});
        std::vector<TrgswDft> trgswDfts(loop, TrgswDft{param});
        Integer mu = 1;
        for (size_t i = 0; i < loop; i++) {
            Integer mui = 3;
            mu *= mui;
            encryptTrgswNtt(trgsws[i], trgswDfts[i], mui, trgswKey, 0, param);
        }
        cout << "mu: " << mu << endl;
        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft{param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu, param.torusBase);
        }
        symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, mu2t);
//        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        Trlwe tmp{param.k, param.N};
        for (size_t i = 0; i < loop; i++) {
            tmp = Trlwe{param.k, param.N};
            externalProductTrgswNtt(tmp, trgswDfts[i], in2, param.lApprox, param);
            swap(in2, tmp);
        }
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, in2, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW_MULT_NAIVE_CHAIN_NTT");
}

TEST(RgswTest, RGSWMP_MULT_NAIVE_CHAIN_NTT) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
    // param.l = 4;
//    param.k = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        int loop = 4;
        std::vector<TrgswMPDft> trgswDfts(loop, TrgswMPDft{param});
        Integer mu = 1;
        for (size_t i = 0; i < loop; i++) {
            Integer mui = 3;
            mu *= mui;
            encryptTrgswMPNtt(trgswDfts[i], mui, trgswKey, 0, param);
        }
        cout << "mu: " << mu << endl;
        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft{param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu, param.torusBase);
        }
        symEncTrlweMultiSampleNtt(in2, in2Dft, trlweKey, mu2t);

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        Trlwe tmp{param.k, param.N};
        for (size_t i = 0; i < loop; i++) {
            tmp = Trlwe{param.k, param.N};
            externalProductTrgswMPNtt(tmp, trgswDfts[i], in2, param.lApprox, param);
            swap(in2, tmp);
        }
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, in2, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_MULT_NAIVE_CHAIN_NTT");
}

// error overflow
TEST(RgswTest, RGSWMP_MULT_NAIVE_DECOMP) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
    param.l = 2;
    // param.lweStdDev = 0;
    // param.rlweStdDev = 0;
//    param.k = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        int loop = 6;
        std::vector<TrgswMP> trgsws(loop, TrgswMP{param});
        Integer mu = 1;
        for (size_t i = 0; i < loop; i++) {
//            Integer mui = 3;
//            mu *= mui;
            encryptLowTrgswMP(trgsws[i], mu, trgswKey, param);
        }
        cout << "mu: " << mu << endl;
        // trlwe enc
        Trlwe in2 {param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = i;
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu, param.torusBase);
        }
        Trlwe out {param.k, param.N};
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
//        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "mu in");

        DecomposedTrlwe outD{param};
        gadgetDecomposeTrlwe(outD, in2, param);

        // trgsw mult
        DecomposedTrlwe tmp{param};
        for (size_t i = 0; i < loop; i++) {
            externalProductTrgswMPDecomp(outD, trgsws[i], outD, param);
        }
//        COUNT_TIME("trgswMPExternalProduct", trgswMPExternalProductDecomp(tmp, trgswMP, in2D, param);)
//        COUNT_TIME("trgswMPExternalProduct", trgswMPExternalProductDecomp(outD, trgswMP2, tmp, param);)
        recomposeTrlwe(out, outD, param);
//        printDecomposedTrlweAB(outD, "outD");
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");
        ASSERT_EQ(multPlain.coeffs, decAftP.coeffs);
    }
    printBanner("RGSWMP_MULT_NAIVE_DECOMP");
}

TEST(RgswTest, RGSW_ROT) {
    YatfheParameters param{};
//    param.N = 1024;
//    param.radixBits = 4;
//    param.l = 3;
//    param.k = 1;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        printMsg(ti, "iter");
        // key gen
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw{param};
//        Integer mu1 = genIntUniformDist(0, 3);
        Integer mu1 = 1;
        int rotN = -512;
        encryptTrgsw(trgsw, mu1, trgswKey, 0, param);
        rotateTrgsw(trgsw, rotN, param);
        printMsg(decryptTrgsw(trgsw, param, trgswKey), "trgsw dec");

        // trlwe enc
        Trlwe in2{param.k, param.N};
        IntPolynomial mu2p{param.N};
        IntPolynomial multPlain{param.N};
        std::vector<Torus> mu2t(param.N);
        for (size_t i = 0; i < param.N; i++) {
            mu2p.coeffs[i] = intModP(i, param.torusBase);
            mu2t[i] = modSwitchToTorus32(mu2p.coeffs[i], param.torusBase);
            multPlain.coeffs[i] = modMulQ(mu2p.coeffs[i], mu1, param.torusBase);
        }
        Trlwe out{param.k, param.N};
        symEncTrlweMultiSample(in2, trlweKey, mu2t);
//        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP{param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult
        COUNT_TIME("trgswExternalProduct", externalProductTrgsw(out, trgsw, in2, param);)
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP{param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");

        // rot back
        TorusPolynomial res{param.N};
        rotateIntPolynomial(res, -rotN, decAftP, param.torusBase);
        printArray(res.coeffs, "res");
        for (auto i = 0 ; i < res.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], res.coeffs[i]);
        }
    }
    printBanner("RGSW_ROT");
}

TEST(RgswTest, RGSW_MULT_MCRT_NAIVE) {
    YatfheParameters param {};
    param.q = Q_CRT;
//    param.N = 32;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        encryptTrgswApproxCRT(trgsw, param, trgswKey, mu1);

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        Integer mu2 = genIntUniformDist(-param.torusBase / 2, (param.torusBase - 1) / 2);
        Torus mu2T = modSwitchToTorus32(mu2, param.torusBase);
        symEncTrlweSingleSample(in2, trlweKey, mu2T);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // GD
        std::vector<Trlwe8> trlweD(param.d, Trlwe8{param.k, param.N});
        decompTrlweMcrt(trlweD, in2, param);
        std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
        decompTrgswMcrt(trgswD, trgsw, param);

        // trgsw mult
        std::vector<Trlwe8> outD(param.d, Trlwe8{param.k, param.N});
        COUNT_TIME("trgswExternalProductApproxCRT", externalProductTrgswApproxCrt(outD, trgswD, trlweD, param);)

        // Recomp
        Trlwe out {param.k, param.N};
        trlweMcrtToCrt(outD, param);
        recompTrlweCrt(out, outD, param);
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        Integer multPlain = intModP(mu1 * mu2, param.torusBase);
        printf("Plain mult: %d * %d = %d\n", mu1, mu2, multPlain);
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain, decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW_MULT_MCRT_NAIVE");
}

TEST(RgswTest, RGSW_MULT_MCRT_NTT) {
    YatfheParameters param {};
    param.q = Q_CRT;
//    param.N = 32;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        encryptTrgswApproxCRT(trgsw, param, trgswKey, mu1);

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        Integer mu2 = genIntUniformDist(-param.torusBase / 2, (param.torusBase - 1) / 2);
        Torus mu2T = modSwitchToTorus32(mu2, param.torusBase);
        symEncTrlweSingleSample(in2, trlweKey, mu2T);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // GD
        std::vector<Trlwe8> trlweD(param.d, Trlwe8{param.k, param.N});
        decompTrlweMcrt(trlweD, in2, param);
        std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
        decompTrgswMcrt(trgswD, trgsw, param);

        // NTT
        std::vector<TrgswDft24> trgswDDft(param.d, TrgswDft24(param.dh, param.k, param.N));
        std::vector<Trgsw8> trgswDI(param.d, Trgsw8(param.dh, param.k, param.N));
        for (size_t d = 0; d < param.d; d++) {
            NttNative24::applyNttForRgsw(trgswDDft[d], trgswD[d]);
        }

        // trgsw mult
        std::vector<Trlwe8> outD(param.d, Trlwe8{param.k, param.N});
        COUNT_TIME("trgswExternalProductApproxCRTNtt", externalProductTrgswApproxCrtNtt(outD, trgswDDft, trlweD, param);)

        // Recomp
        Trlwe out {param.k, param.N};
        trlweMcrtToCrt(outD, param);
        recompTrlweCrt(out, outD, param);
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        Integer multPlain = intModP(mu1 * mu2, param.torusBase);
        printf("Plain mult: %d * %d = %d\n", mu1, mu2, multPlain);
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain, decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW_MULT_MCRT_NTT");
}

//todo: trlgsw
TEST(RgswTest, RgswMultTestNTT14) {
    YatfheParameters param {};
    param.radixBits = 8;
//    param.l = 4;
    param.l2 = 4;
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
//        Integer mu1 = genIntUniformDist(1, 3);
        Integer mu1 = 3;
        encryptTrgswNtt(trgsw, trgswDft, mu1, trgswKey, 0, param);

        Trlgsw trlgsw {param};
        TrlgswDft14 trlgswDft14 {param};
        encryptTrlgswNtt(trlgsw, trlgswDft14, param, trgswKey, mu1);

        printf( "trgsw dec: %d.\n", decryptTrgswNtt(trgswDft, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft {param.k, param.N};
//        Integer mu2p = genIntUniformDist(1, 1);
        Integer mu2p = 1;
        Torus mu2 = modSwitchToTorus32(mu2p, param.torusBase);
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        symEncTrlweSingleSampleNtt(in2, in2Dft, trlweKey, mu2);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        // todo
        Trlwe out {param.k, param.N};
        Trlwe out14p {param.k, param.N};
        Trlwe out14 {param.k, param.N};
        COUNT_TIME("trgswExternalProductNtt", externalProductTrgswNtt(out, trgswDft, in2, param.lApprox, param);)
        externalProductTrlgsw(out14p, trlgsw, in2, param);
        COUNT_TIME("trlgswExternalProductNtt14", externalProductTrlgswNtt(out14, trlgswDft14, in2, param);)

        printTrlweAB(out, "out");
        printTrlweAB(out14p, "out14p");
        printTrlweAB(out14, "out14");


        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(intModP(mu1 * mu2p, param.torusBase), decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNTT");
}