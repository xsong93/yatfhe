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

    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    lweKeyGen(tlweKey);
    trlweKeyGen(trlweKey);

    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Integer plain  = 7;
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
    Torus dec = trgswDecrypt(trgswDft, param, trgswKey);
    cout << "plain: " << plain << endl;
    cout << "dec: " << dec << endl;
    ASSERT_EQ(plain, dec);
}

TEST(RgswMultTest, RgswMultTest) {
    const YatfheParameters param {};

    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    lweKeyGen(tlweKey);
    trlweKeyGen(trlweKey);

    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    Torus mu1 = doubleToTorus32(1.0 / param.torusBase);
    trgswEncrypt(trgsw, trgswDft, param, trgswKey, mu1);
//    printTrgsw(trgsw, "trgsw");

    Trlwe in2 {param.k, param.N};
    TrlweDft in2Dft {param.k, param.N};
    Torus mu2 = doubleToTorus32(1.0 / param.torusBase);
    Trlwe out {param.k, param.N};
    symEncTrlweSingleSample(in2, in2Dft, trlweKey, mu2, param.lweStdDev);
    // todo: mult
}