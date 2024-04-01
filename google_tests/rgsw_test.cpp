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
    lweKeyGen(tlweKey, param.n);
    setCoeffsValue(tlweKey.s, 0); // todo: remove this
    trlweKeyGen(trlweKey, param.N, param.k);
    for (auto i = 0; i < trlweKey.k; i++) {
        setCoeffsValue(trlweKey.s[i].coeffs, 0); // todo: remove this
    }
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    trgswEncrypt(trgsw, trgswDft, param, trgswKey, 0);
    printTrgsw(trgsw, "trgsw");

    // test identity for trgsw and trgswDft value
    for (auto i = 0; i < trgsw.l; i++) {
        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
            TorusPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
            ASSERT_EQ(ip.coeffs, trgsw.trlweSamples[i][j].b.coeffs);
        }
    }
    //todo: test on cmux
    Trlwe in {param.k, param.N};
    Trlwe out {param.k, param.N};
    for (auto j = 0; j < in.b.N; j++) {
        for (auto i = 0; i < in.k; i++) {
            in.a[i].coeffs[j] = j;
        }
        in.b.coeffs[j] = j;
    }
    controlMux(out, in, 1, trgswDft, param);
    printTrlweAB(in, "in");
    printTrlweAB(out, "out");

}

TEST(RgswMultTest, RgswMultTest) {
    const YatfheParameters param {};

    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    lweKeyGen(tlweKey, param.n);
    trlweKeyGen(trlweKey, param.N, param.k);

    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    trgswEncrypt(trgsw, trgswDft, param, trgswKey, doubleToTorus32(1.0 / param.torusBase));
//    printTrgsw(trgsw, "trgsw");

    // test identity for trgsw and trgswDft value
    for (auto i = 0; i < trgsw.l; i++) {
        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
            TorusPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
            ASSERT_EQ(ip.coeffs, trgsw.trlweSamples[i][j].b.coeffs);
        }
    }
    Trlwe in2 {param.k, param.N};
    Trlwe out {param.k, param.N};
    // todo: mult
}