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
    setCoeffsValue(tlweKey.s, 1); // todo: remove this
    trlweKeyGen(trlweKey, param.N, param.k);
    for (auto i = 0; i < trlweKey.k; i++) {
        setCoeffsValue(trlweKey.s[i].coeffs, 1); // todo: remove this
    }
    Trgsw trgsw {param};
    TrgswDft trgswDft {param};
    trgswEncrypt(trgsw, trgswDft, param, trgswKey, 1);
    for (auto i = 0; i < trgsw.l; i++) {
        for (auto j = 0; j < trgsw.trlweSamples[i].size(); j++) {
            IntPolynomial ip {trgswDft.trlweDftSamples[i][j].b.N};
            applyIntt(ip, trgswDft.trlweDftSamples[i][j].b);
//            printTrlweAB(trgsw.trlweSamples[i][j], "trgsw.trlweSamples");
            printArray(trgsw.trlweSamples[i][j].b.coeffs, "ori");
            printArray(ip.coeffs, "ntt");
        }
    }
}