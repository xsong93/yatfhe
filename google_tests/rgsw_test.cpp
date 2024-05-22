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

TEST(RgswTest, RgswEncDecTest) {
    const YatfheParameters param {};
    initGlobalParamsNtt64(param.N);

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
//    Integer dec = trgswDecryptNtt(trgswDft, param, trgswKey);
    Integer dec = trgswDecrypt(trgsw, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
    printBanner("RgswEncDecTest");
}

TEST(RgswMultTest, RgswMultTest) {
    YatfheParameters param {};
    param.N = 1024;
    param.radixBits = 4;
    param.l = 8;
    param.k = 2;
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
        Integer mu1 = genIntUniformDist(0, 1);
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
        symDecTrlweToInt(decAftP, out, trlweKey, param.torusBase);
        printArray(decAftP.coeffs, "decAftP");
        for (auto i = 0 ; i < decAftP.N; i++) {
            ASSERT_EQ(intModP(mu1 * mu2p, param.torusBase), decAftP.coeffs[i]);
        }
    }
    printBanner("RgswMultTest");
}