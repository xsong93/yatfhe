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
#include "yatfhe/ntt24.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

TEST(RgswTest, RgswEncDecTest) {
    YatfheParameters param {};
    yatfheInit(param);

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain = 7;
    trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, plain);
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
    Integer dec = trgswDecryptNtt(trgswDft, param, trgswKey);
//    Integer dec = trgswDecrypt(trgsw, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RgswEncDecTest");
}

TEST(RgswTest, RGSW_MCRT_DECOMPOSITION) {
    YatfheParameters param {};
    yatfheInit(param);

    // key gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain = 7;
    trgswEncrypt(trgsw, param, trgswKey, plain);

    // GD
    std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
    Trgsw trgswRe {param};
    trgswMCRTDecomp(trgswD, trgsw, param);
    trgswMCRTToCRT(trgswD, param);
    trgswCRTRecomp(trgswRe, trgswD, param);
    printTrgsw(trgsw, "trgsw");
    printTrgsw(trgswRe, "trgswRe");

    // trgsw dec
    Integer dec = trgswDecrypt(trgswRe, param, trgswKey);
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
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

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
        COUNT_TIME("trgswExternalProduct", trgswExternalProduct(out, trgsw, in2, param);)
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

TEST(RgswTest, RGSWMP_MULT_NAIVE) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
//    param.l = 2;
//    param.k = 3;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        TrgswMP trgswMP {param};
        Integer mu1 = 1;
        trgswMPEncrypt(trgswMP, mu1, param, trgswKey);

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
        COUNT_TIME("trgswMPExternalProduct", trgswMPExternalProduct(out, trgswMP, in2, param);)
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

// chaining multiplication only works for exact decomp when multiplier is greater than 1
TEST(RgswTest, RGSWMP_MULT_NAIVE_CHAIN) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
    param.l = 4;
//    param.k = 3;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        int loop = param.n;
        std::vector<Trgsw> trgsws(loop, Trgsw{param});
        std::vector<TrgswDft> trgswDfts(loop, TrgswDft{param});
        Integer mu = 1;
        for (size_t i = 0; i < loop; i++) {
            Integer mui = 3;
            mu *= mui;
            trgswEncryptNtt(trgsws[i], trgswDfts[i], param, trgswKey, mui);
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
            trgswExternalProductNtt(tmp, trgswDfts[i], in2, param);
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
    printBanner("RGSWMP_MULT_NAIVE_CHAIN");
}

// chaining multiplication only works for exact decomp when multiplier is greater than 1
TEST(RgswTest, RGSWMP_MULT_NAIVE_DECOMP) {
    YatfheParameters param {};
//    param.n = 805;
//    param.N = 512;
//    param.radixBits = 10;
    param.l = 4;
//    param.k = 3;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        int loop = param.n;
        std::vector<TrgswMP> trgsws(loop, TrgswMP{param});
        Integer mu = 1;
        for (size_t i = 0; i < loop; i++) {
            Integer mui = 3;
            mu *= mui;
            trgswMPEncryptLow(trgsws[i], mui, param, trgswKey);
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
            tmp = DecomposedTrlwe{param};
            trgswMPExternalProductDecomp(tmp, trgsws[i], outD, param);
            swap(outD, tmp);
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
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], decAftP.coeffs[i]);
        }
    }
    printBanner("RGSWMP_MULT_NAIVE_DECOMP");
}

TEST(RgswTest, RGSW_ROT) {
    YatfheParameters param{};
//    param.N = 1024;
//    param.radixBits = 4;
//    param.l = 3;
//    param.k = 1;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw{param};
//        Integer mu1 = genIntUniformDist(0, 3);
        Integer mu1 = 1;
        int rotN = 513;
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        COUNT_TIME("trgswRotate", trgswRotate(trgsw, rotN, param);)
        printf("trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

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
        COUNT_TIME("trgswExternalProduct", trgswExternalProduct(out, trgsw, in2, param);)
//        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        IntPolynomial decAftP{param.N};
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        printArray(multPlain.coeffs, "Plain mult");

        // rot back
        TorusPolynomial res{param.N};
        intPolynomialRotate(res, -rotN, decAftP, param.torusBase);
        printArray(res.coeffs, "res");
        for (auto i = 0 ; i < res.N; i++) {
            ASSERT_EQ(multPlain.coeffs[i], res.coeffs[i]);
        }
    }
    printBanner("RGSW_ROT");
}

TEST(RgswTest, RGSW_MULT_MCRT_NAIVE) {
    YatfheParameters param {};
//    param.N = 32;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        trgswEncryptApproxCRT(trgsw, param, trgswKey, mu1);

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
        trlweMCRTDecomp(trlweD, in2, param);
        std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
        trgswMCRTDecomp(trgswD, trgsw, param);

        // trgsw mult
        std::vector<Trlwe8> outD(param.d, Trlwe8{param.k, param.N});
        COUNT_TIME("trgswExternalProductApproxCRT", trgswExternalProductApproxCRT(outD, trgswD, trlweD, param);)

        // Recomp
        Trlwe out {param.k, param.N};
        trlweMCRTToCRT(outD, param);
        trlweCRTRecomp(out, outD, param);
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
//    param.N = 32;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        trgswEncryptApproxCRT(trgsw, param, trgswKey, mu1);

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
        trlweMCRTDecomp(trlweD, in2, param);
        std::vector<Trgsw8> trgswD(param.d, Trgsw8(param.dh, param.k, param.N));
        trgswMCRTDecomp(trgswD, trgsw, param);

        // NTT
        std::vector<TrgswDft24> trgswDDft(param.d, TrgswDft24(param.dh, param.k, param.N));
        std::vector<Trgsw8> trgswDI(param.d, Trgsw8(param.dh, param.k, param.N));
        for (size_t d = 0; d < param.d; d++) {
            applyNttForRgsw24(trgswDDft[d], trgswD[d]);
        }

        // trgsw mult
        std::vector<Trlwe8> outD(param.d, Trlwe8{param.k, param.N});
        COUNT_TIME("trgswExternalProductApproxCRTNtt", trgswExternalProductApproxCRTNtt(outD, trgswDDft, trlweD, param);)

        // Recomp
        Trlwe out {param.k, param.N};
        trlweMCRTToCRT(outD, param);
        trlweCRTRecomp(out, outD, param);
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

TEST(RgswTest, RgswMultTestNTT) {
    YatfheParameters param {};
//    param.k = 8;
//    param.N = 128;
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
        Integer mu1 = genIntUniformDist(0, 3);
        trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecryptNtt(trgswDft, param, trgswKey));

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
        COUNT_TIME("trgswExternalProductNtt", trgswExternalProductNtt(out, trgswDft, in2, param);)
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

//todo: trlgsw
TEST(RgswTest, RgswMultTestNTT14) {
    YatfheParameters param {};
    param.radixBits = 8;
    param.l = 4;
    param.l2 = 4;
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
//        Integer mu1 = genIntUniformDist(1, 3);
        Integer mu1 = 3;
        trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, mu1);

        Trlgsw trlgsw {param};
        TrlgswDft14 trlgswDft14 {param};
        trlgswEncryptNtt14(trlgsw, trlgswDft14, param, trgswKey, mu1);

        printf( "trgsw dec: %d.\n", trgswDecryptNtt(trgswDft, param, trgswKey));

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
        COUNT_TIME("trgswExternalProductNtt", trgswExternalProductNtt(out, trgswDft, in2, param);)
        trlgswExternalProduct(out14p, trlgsw, in2, param);
        COUNT_TIME("trlgswExternalProductNtt14", trlgswExternalProductNtt14(out14, trlgswDft14, in2, param);)

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