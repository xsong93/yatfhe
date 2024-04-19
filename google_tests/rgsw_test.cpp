//
// Created by Xintong Song on 2024/3/26.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yatfhe/ntt.h"

TEST(RgswEncDecTest, RgswEncDecTest) {
    const YatfheParameters param {};

    // ken gen
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    trlweKeyGen(trlweKey);

    // trgsw enc
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain = 7;
//    trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, plain);
    trgswEncrypt(trgsw, param, trgswKey, plain);
//    printTrgsw(trgsw, "trgsw");

//    // test identity for trgsw and trgswDft value
//    for (auto i = 0; i < trgsw.l; i++) {
//        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
//            TorusPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
//            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
//            ASSERT_EQ(ip.coeffs, trgsw.trlweSamples[i][j].b.coeffs);
//        }
//    }

    // trgsw dec
//    Torus dec = trgswDecryptNtt(trgswDft, param, trgswKey);
    Torus dec = trgswDecrypt(trgsw, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
}

TEST(RgswMultTest, RgswMultTest) {
    YatfheParameters param {};
    param.N = 1;
    param.radixBits = 4;
    param.l = 8;
    param.k = 2;
    int ti = 1;
    while (ti-- > 0) {
        // ken gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
        Integer mu1 = 1;
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        Torus mu2 = doubleToTorus32(1.0 / param.torusBase);
        Trlwe out {param.k, param.N};
        DoublePolynomial decPre {param.N};
        DoublePolynomial decAft {param.N};
        symEncTrlweSingleSample(in2, trlweKey, mu2, param.rlweStdDev);
        printTrlweAB(in2, "trlwe");

        // trlwe dec pre-mult
        symDecTrlwe(decPre, in2, trlweKey, param.torusBase);
        printArray(decPre.coeffs, "decPre");

        // trgsw mult
        // todo: ntt/intt failed, fix ntt
        trgswExternalProduct(out, trgsw, in2, param);
        printTrlweAB(out, "out");

//        // intt and test
//        applyNttForAB(in2Dft, out);
//        Trlwe out2 {param.k, param.N};
//        applyInttForAB(out2, in2Dft);
//        for (auto i = 0; i < out2.b.N; i++) {
//            for (auto j = 0; j < out2.k; j++) {
//                ASSERT_EQ(out2.a[j].coeffs[i], out.a[j].coeffs[i]);
//            }
//            ASSERT_EQ(out2.b.coeffs[i], out.b.coeffs[i]);
//        }

        // trlwe dec aft-mult
        symDecTrlwe(decAft, out, trlweKey, param.torusBase);
        printArray(decAft.coeffs, "decAft");
        for (auto i = 0 ; i < decAft.N; i++) {
            ASSERT_EQ(decPre.coeffs[i], decAft.coeffs[i]);
        }
    }
}