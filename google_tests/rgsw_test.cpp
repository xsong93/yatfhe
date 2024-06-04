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

// todo
TEST(RgswTest, RgswEncDecTestNtt14) {
    YatfheParameters param {};
    yatfheInit(param);

    // ken gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    TrlgswDft14 trlgswDft14 {param};
    Integer plain = 7;

    trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, plain);
    trgswEncryptNtt14(trgsw, trlgswDft14, param, trgswKey, plain);
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
    printBanner("RgswEncDecTestNtt14");
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
        TrgswDft trgswDft {param};
        Integer mu1 = genIntUniformDist(0, 3);
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        Integer mu2p = genIntUniformDist(INT32_MIN, INT32_MAX);
        Torus mu2 = modSwitchToTorus32(mu2p, param.torusBase);
        Trlwe out {param.k, param.N};
//        DoublePolynomial decPre {param.N};
//        DoublePolynomial decAft {param.N};
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        symEncTrlweSingleSample(in2, trlweKey, mu2);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlweToInt(decPreP, in2, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult
        COUNT_TIME("trgswExternalProduct", trgswExternalProduct(out, trgsw, in2, param);)
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(intModP(mu1 * mu2p, param.torusBase), decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNaive");
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
        Integer mu2p = genIntUniformDist(INT32_MIN, INT32_MAX);
        Torus mu2 = modSwitchToTorus32(mu2p, param.torusBase);
        Trlwe out {param.k, param.N};
        IntPolynomial decPreP {param.N};
        IntPolynomial decAftP {param.N};
        symEncTrlweSingleSampleNtt(in2, in2Dft, trlweKey, mu2);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlweToIntNtt(decPreP, in2Dft, trlweKey, param.torusBase);
        printArray(decPreP.coeffs, "decPreP");

        // trgsw mult ntt
        // todo
        COUNT_TIME("trgswExternalProductNtt", trgswExternalProductNtt(out, trgswDft, in2, param);)
        printTrlweAB(out, "out");

        // trlwe dec aft-mult
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(intModP(mu1 * mu2p, param.torusBase), decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTestNTT");
}