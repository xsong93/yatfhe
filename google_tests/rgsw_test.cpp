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
    const YatfheParameters param {};

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
    DoublePolynomial dec {param.N};
    symEncTrlweSingleSample(in2, in2Dft, trlweKey, mu2, param.lweStdDev);
    printTrlweAB(in2, "in2");

    // trlwe dec pre-mult
    symDecTrlwe(dec, in2Dft, trlweKey, param.torusBase);
    printArray(dec.coeffs, "decPre");

    // trgsw mult
    trgswExternalProduct(out, trgswDft, in2, param);
    applyNttForAB(in2Dft, out);

    // trlwe dec aft-mult
    symDecTrlwe(dec, in2Dft, trlweKey, param.torusBase);
    printArray(dec.coeffs, "decAft");
}