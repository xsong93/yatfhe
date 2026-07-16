//
// Created by Xintong Song on 2024/3/26.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlgsw.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/ntt_hexl.h"

TEST(RGSW, ENC_DEC) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    auto level = param.lApprox;

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw {param, level};
    TrgswDft trgswDft {param, level};
    Integer plain = 1;
    encryptTrgswNtt(trgsw, trgswDft, plain, trgswKey, 0, param);

    // trgsw dec
    Integer dec = decryptTrgswNtt(trgswDft, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RGSW.ENC_DEC");
}

TEST(RGSW, ADD) {
    YatfheParameters param{};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    Trgsw trgsw1{param}, trgsw2{param}, res{param};
    TrgswDft trgswDft1{param}, trgswDft2{param}, resDft{param};
    Integer dec, plainAdd;
    int plain1, plain2;
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
    printBanner("RGSW.ADD");
}

TEST(RGSW, SUB) {
    YatfheParameters param{};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
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
    printBanner("RGSW.SUB");
}

TEST(RGSW, MCRT_DECOMPOSITION) {
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
    // printTrgsw(trgsw, "trgsw");
    // printTrgsw(trgswRe, "trgswRe");

    // trgsw dec
    Integer dec = decryptTrgsw(trgswRe, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RGSW.MCRT_DECOMPOSITION");
}

TEST(RGSW, MULT_NAIVE) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
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
        // printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "decPreP");

        // trgsw mult
        COUNT_TIME("trgswExternalProduct", externalProductTrgsw(out, trgsw, in2, param);)
        // printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");
        // printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW.MULT_NAIVE");
}

TEST(RGSW, MULT_NTT) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
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
        // printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        COUNT_TIME("trgswExternalProductNtt", externalProductTrgswNtt(out, trgswDft, in2, param.lApprox, param);)
        // printTrlweAB(out, "out");

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");

        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        // printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW.MULT_NTT");
}

TEST(RGSWMP, ENC_DEC) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    auto level = param.lApprox;

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    TrgswMPDft trgswDft {param, level};
    int plain = 1;
    int pos = 0;
    encryptTrgswMPNtt(trgswDft, plain, trgswKey, pos, param);

    // trgsw dec
    IntPolynomial dec {param.N};
    decryptTrgswMPNtt(dec, trgswDft, param, trgswKey, false);

    for (auto i = 0; i < dec.N; i++) {
        if (i == pos) ASSERT_EQ(plain, dec.coeffs[i]);
        else ASSERT_EQ(0, dec.coeffs[i]);
    }
    printBanner("RGSWMP.ENC_DEC");
}

TEST(RGSWMP, MULT_NAIVE) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc
        TrgswMP trgswMP {param};
        Integer mu1 = genIntUniformDist(0, 3);
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

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "mu in");

        // trgsw mult
        COUNT_TIME("trgswMPExternalProduct", externalProductTrgswMP(out, trgswMP, in2, param.lApprox, param);)

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");
        // printArray(multPlain.coeffs, "Plain mult");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP.MULT_NAIVE");
}

TEST(RGSWMP, MULT_CONST) {
    YatfheParameters param {};
    param.torusBase = 16;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    int ti = 0;
    while (ti++ < 8) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc

        Integer mu1 = genIntUniformDist(-4, 3);

        TrgswMP trgswMP {param};
        vector<Integer> mu2ps(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2ps[i] = genIntUniformDist(-1, 1);
        }
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        encryptTrgswMPMulti(trgswMP, mu2ps, trgswKey, param);

        // trgsw dec pre-mult
        decryptTrgswMP(decPreP, trgswMP, param, trgswKey, false);
        // printArray(decPreP.coeffs, "decPreP");

        multTrgswMPWithConst(trgswMP, trgswMP, mu1);

        decryptTrgswMP(decAftP, trgswMP, param, trgswKey, false);
        // printArray(decAftP.coeffs, "decAftP");


        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        // printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP.MULT_CONST");
}

TEST(RGSWMP, MULT_CONST_NTT) {
    YatfheParameters param {};
    param.torusBase = 16;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    int ti = 0;
    while (ti++ < 8) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTrlweKey(trlweKey);

        // trgsw enc

        Integer mu1 = genIntUniformDist(-4, 3);

        TrgswMPDft trgswMP {param};
        vector<Integer> mu2ps(param.N);
        for (auto i = 0; i < param.N; i++) {
            mu2ps[i] = genIntUniformDist(-1, 1);
        }
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        encryptTrgswMPMultiNtt(trgswMP, mu2ps, trgswKey, param);

        // trgsw dec pre-mult
        decryptTrgswMPNtt(decPreP, trgswMP, param, trgswKey, false);
        // printArray(decPreP.coeffs, "decPreP");

        multTrgswMPWithConstNtt(trgswMP, trgswMP, mu1);

        decryptTrgswMPNtt(decAftP, trgswMP, param, trgswKey, false);
        // printArray(decAftP.coeffs, "decAftP");


        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        // printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP.MULT_CONST_NTT");
}

TEST(RGSWMP, MULT_NTT) {
    YatfheParameters param {};
    param.torusBase = 16;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
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
        Trlwe out {param};

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        COUNT_TIME("trgswMPExternalProductNtt", externalProductTrgswMPNtt(out, trgswDft, in2, param.lApprox, param);)

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");

        vector<Integer> pMult(param.N);
        for (auto i = 0; i < param.N; i++) {
            pMult[i] = intModP(mu1 * mu2ps[i], param.torusBase);
        }
        // printArray(pMult, "realVal");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(pMult[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP.MULT_NTT");
}

TEST(RGSWMP, INTERMULT_NAIVE) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    int level = param.l;

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    TrgswMP trgswMP1 {param, level};
    Integer mu1 = genIntUniformDist(0, 3);
    encryptTrgswMP(trgswMP1, mu1, trgswKey, 0, param);

    TrgswMP trgswMP2 {param, level};
    Integer mu2 = genIntUniformDist(0, 3);
    encryptTrgswMP(trgswMP2, mu2, trgswKey, 0, param);

    // trgsw mult
    Trlwe out{param.k, param.N};
    IntPolynomial t1{param.N};
    // decryptTrgswMP(t1, trgswMP1, param, trgswKey, false);
    // printArray(t1.coeffs, "trgswMP1");
    // decryptTrgswMP(t1, trgswMP2, param, trgswKey, false);
    // printArray(t1.coeffs, "trgswMP2");

    TrgswMP mult{param, level};
    COUNT_TIME("trgswMPInternalProduct",
        internalProductTrgswMP(mult, trgswMP1, trgswMP2, level, param);)

    decryptTrgswMP(t1, mult, param, trgswKey, false);
    // printArray(t1.coeffs, "dec");
    ASSERT_EQ(intModP(mu1*mu2, param.torusBase), t1.coeffs[0]);
    for (auto i = 1 ; i < t1.N; i++) {
        ASSERT_EQ(0, t1.coeffs[i]);
    }

    printBanner("RGSWMP.INTERMULT_NAIVE");
}

TEST(RGSWMP, INTERMULT_NAIVE_NTT) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    TrgswMP trgswMP1{param};
    Integer mu1 = genIntUniformDist(0, 3);
    encryptTrgswMP(trgswMP1, mu1, trgswKey, 0, param);

    TrgswMPDft trgswMP2{param};
    Integer mu2 = genIntUniformDist(0, 3);
    encryptTrgswMPNtt(trgswMP2, mu2, trgswKey, 0, param);

    TrgswMP mult{param};
    COUNT_TIME("internalProductTrgswMPNtt",
        internalProductTrgswMPNtt(mult, trgswMP1, trgswMP2, param.l, param);)

    IntPolynomial dec{param.N};
    decryptTrgswMP(dec, mult, param, trgswKey, false);
    // printArray(dec.coeffs, "dec");
    ASSERT_EQ(intModP(mu1*mu2, param.torusBase), dec.coeffs[0]);
    for (auto i = 1 ; i < dec.N; i++) {
        ASSERT_EQ(0, dec.coeffs[i]);
    }
    printBanner("RGSWMP.INTERMULT_NAIVE_NTT");
}

TEST(RGSWMP, SCHEME_SWITCHING) {
    YatfheParameters param{};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey &trlweKey = trgswKey.trlweKey;
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
    {
        // encryptTrgswMPNtt samples cPrime's "a" natively/uniformly over the wide qNtt domain,
        // which switchTrlweToSecretEmbeddingNtt cannot correctly INTT+decompose (it needs a
        // genuinely bounded Torus-domain "a" forward-transformed into NTT domain instead).
        TorusPolynomial muPolyFix{param.N};
        Trlwe scratch{param.k, param.N};
        for (auto lvl = 0; lvl < in1Dft.l; lvl++) {
            muPolyFix.coeffs[0] = static_cast<Torus>(mu1) << (param.torusBits - (lvl + 1) * param.radixBits);
            symEncTrlweMultiSampleNtt(scratch, in1Dft.cPrime[lvl], trgswKey.trlweKey, muPolyFix.coeffs);
        }
    }

    IntPolynomial dec{param.N};
    decryptTrgswMPNtt(dec, in1Dft, param, trgswKey, false);
    // printArray(dec.coeffs, "mu1");

    // init data
    TrgswMPDft tmpMp{param, param.lApprox};
    tmpMp.cPrime = in1Dft.cPrime;
    TrgswMPDft tmpMp2{param, param.lApprox};
    vector decompA(param.lApprox, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    for (auto l = 0; l < param.lApprox; l++) {
        NttHexl::applyNtt(tmpMp2.cPrime[l].b, in1.cPrime[l].b);
    }
    for (auto l0 = 0; l0 < param.lApprox; l0++) {
        auto &a = in1.cPrime[l0].a;
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
    TrgswMPDft tmpMp3{param, param.lApprox};
    TrgswMP t1{param, param.lApprox};
    t1.cPrime = in1.cPrime;

    // scheme switching
    for (auto l = 0; l < param.lApprox; l++) {
        COUNT_TIME("switchTrlweToSecretEmbedding", switchTrlweToSecretEmbedding(t1.c[l], in1.cPrime[l], s2p, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNtt",
                   switchTrlweToSecretEmbeddingNtt(tmpMp.c[l], tmpMp.cPrime[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNttOpt",
                   switchTrlweToSecretEmbeddingNttOpt(tmpMp2.c[l], tmpMp2.cPrime[l], decompA[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNttMix",
                   switchTrlweToSecretEmbeddingNttMix(tmpMp3.c[l], tmpMp3.cPrime[l], decompA[l], in1.cPrime[l].b,
                                                      s2pDft, param);)
    }

    // result validation
    Trlwe trlwe{param.k, param.N};
    Trlwe res{param.k, param.N};
    TorusPolynomial mu{param.N};
    IntPolynomial plainMult{param.N};
    symEncTrlweSingleSample(trlwe, trlweKey, modSwitchToTorus32(3, param.torusBase), 0);

    // switchTrlweToSecretEmbedding
    externalProductTrgswMP(res, t1, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    // printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbedding: " << (dec.coeffs[0] == 3) << endl;

    // switchTrlweToSecretEmbeddingNtt
    COUNT_TIME("externalProductTrgswMPNtt", externalProductTrgswMPNtt(res, tmpMp, trlwe, param.lApprox, param);)
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    // printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNtt: "  << (dec.coeffs[0] == 3) << endl;

    // switchTrlweToSecretEmbeddingNttOpt
    externalProductTrgswMPNtt(res, tmpMp2, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    // printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNttOpt: "  << (dec.coeffs[0] == 3) << endl;

    // switchTrlweToSecretEmbeddingNttMix
    externalProductTrgswMPNtt(res, tmpMp3, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    // printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNttMix: "  << (dec.coeffs[0] == 3) << endl;

    printBanner("RGSWMP.SCHEME_SWITCHING");
}

TEST(RGSW, ROT) {
    YatfheParameters param{};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
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
        // printMsg(decryptTrgsw(trgsw, param, trgswKey), "trgsw dec");

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
        // printArray(decPreP.coeffs, "decPreP");

        // trgsw mult
        COUNT_TIME("trgswExternalProduct", externalProductTrgsw(out, trgsw, in2, param);)
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP{param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");
        // printArray(multPlain.coeffs, "Plain mult");

        // rot back
        TorusPolynomial res{param.N};
        rotateIntPolynomial(res, -rotN, decAftP, param.torusBase);
        // printArray(res.coeffs, "res");
        for (auto i = 0 ; i < res.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], res.coeffs[i]);
        }
    }
    printBanner("RGSW.ROT");
}

TEST(RGSW, MULT_MCRT_NAIVE) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.q = Q_CRT;
    param.l = param.d;
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
        symEncTrlweSingleSample(in2, trlweKey, mu2T, 0);
        // printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "decPreP");

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
        // printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");
        Integer multPlain = intModP(mu1 * mu2, param.torusBase);
        // printf("Plain mult: %d * %d = %d\n", mu1, mu2, multPlain);
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(i == 0 ? multPlain : 0, decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW.MULT_MCRT_NAIVE");
}

TEST(RGSW, MULT_MCRT_NTT) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.q = Q_CRT;
    param.l = param.d;
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
        symEncTrlweSingleSample(in2, trlweKey, mu2T, 0);
        // printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        IntPolynomial decPreP {param.N};
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        // printArray(decPreP.coeffs, "decPreP");

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
        // printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP {param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        // printArray(decAftP.coeffs, "decAftP");
        Integer multPlain = intModP(mu1 * mu2, param.torusBase);
        // printf("Plain mult: %d * %d = %d\n", mu1, mu2, multPlain);
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(i == 0 ? multPlain : 0, decAftP.coeffs[i]);
        }
    }
    printBanner("RGSW.MULT_MCRT_NTT");
}