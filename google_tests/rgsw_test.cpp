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
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

TEST(RgswTest, RgswEncDecTest) {
    YatfheParameters param {};
    yatfheInit(param);

    // ken gen
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

    // ken gen
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
        // ken gen
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
        Integer mu2p = genIntUniformDist(-param.torusBase / 2, (param.torusBase - 1) / 2);
        Torus mu2 = modSwitchToTorus32(mu2p, param.torusBase);
        Trlwe out {param.k, param.N};
        symEncTrlweSingleSample(in2, trlweKey, mu2);
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
        Integer multPlain = intModP(mu1 * mu2p, param.torusBase);
        printf("Plain mult: %d * %d = %d\n", mu1, mu2p, multPlain);
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(multPlain, decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNaive");
}

//todo
TEST(RgswTest, RGSW_MULT_MCRT_NAIVE) {
    YatfheParameters param {};
//    param.N = 1024;
//    param.radixBits = 4;
//    param.l = 3;
//    param.k = 1;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 1) {
        cout << "iter: " << ti << endl;
        // ken gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        Integer mu1 = genIntUniformDist(0, 3);
        mu1 = 1;
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

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
        COUNT_TIME("trgswExternalProductCRT", trgswExternalProductCRT(outD, trgswD, trlweD, param);)

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

TEST(RgswTest, RgswMultTestNTT) {
    YatfheParameters param {};
    printf("n:%d, k:%d, N:%d, b:%d, l:%d", param.n, param.k, param.N, param.radixBits, param.l);
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;
        // ken gen
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
        // ken gen
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