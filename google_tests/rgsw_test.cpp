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
    trgswEncrypt(trgsw, trgswDft, param, trgswKey, plain);
//    printTrgsw(trgsw, "trgsw");

    // test identity for trgsw and trgswDft value
    for (auto i = 0; i < trgsw.l; i++) {
        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
            TorusPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
            ASSERT_EQ(ip.coeffs, trgsw.trlweSamples[i][j].b.coeffs);
        }
    }

    // trgsw dec
    Torus dec = trgswDecrypt(trgswDft, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
}

TEST(RgswMultTest, RgswMultTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.l = 8;
    param.k = 2;
    int ti = 20;
    while (ti-- > 0) {
        // ken gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
        Integer mu1 = 1;
        trgswEncrypt(trgsw, trgswDft, param, trgswKey, mu1);
        printf( "trgsw: %d.\n", trgswDecrypt(trgswDft, param, trgswKey));

        // trlwe enc
        Trlwe in2 {param.k, param.N};
        TrlweDft in2Dft {param.k, param.N};
        Torus mu2 = doubleToTorus32(1.0 / param.torusBase);
        Trlwe out {param.k, param.N};
        DoublePolynomial decPre {param.N};
        DoublePolynomial decAft {param.N};
        symEncTrlweSingleSample(in2, in2Dft, trlweKey, mu2, param.lweStdDev);
//        printTrlweAB(in2, "in2");

        // trlwe dec pre-mult
        symDecTrlwe(decPre, in2Dft, trlweKey, param.torusBase);
//        printArray(decPre.coeffs, "decPre");

        // trgsw mult
        // todo: ntt/intt failed, fix ntt
        trgswExternalProduct(out, trgswDft, in2, param);

        Trlwe minus {param.k, param.N};
        trlweSub(minus, in2, out);
        printTrlweAB(minus, "minus");
        for (auto i = 0 ; i < out.b.N; i++) {
//            for (auto j = 0 ; j < out.k; j++) {
//                if (minus.a[j].coeffs[i] != 0) {
//                    printf("%d,%d\n", j, i);
//                }
//                ASSERT_EQ(minus.a[j].coeffs[i], 0);
//            }
            if (minus.b.coeffs[i] != 0) {
                printf("%d\n",i);
            }
            ASSERT_EQ(minus.b.coeffs[i], 0);
        }

//        applyNttForAB(in2Dft, out);
//
//        // trlwe dec aft-mult
//        symDecTrlwe(decAft, in2Dft, trlweKey, param.torusBase);
//        printArray(decAft.coeffs, "decAft");
//        for (auto i = 0 ; i < decAft.N; i++) {
//            ASSERT_EQ(decPre.coeffs[i], decAft.coeffs[i]);
//        }
    }
}