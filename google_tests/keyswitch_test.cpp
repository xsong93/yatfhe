//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(KSKTest, KSKTest) {
    YatfheParameters param {};
    param.torusBase = 128;
    initYatfhe(param);
    TlweKey tlweKey {param.n, param.lweStdDev};
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TlweKey inKey(param.k * param.N, param.lweStdDev);
    convertTrlweKeyToTlweKey(inKey, trlweKey);
    for (int i = - param.torusBase / 2; i < param.torusBase / 2; i++) {
        Integer plainMsg = i;
        Torus mu = modSwitchToTorus32(plainMsg, param.torusBase);
        Tlwe input {param.N * param.k};
        Tlwe output {param.n};
        symEncTlwe(input, mu, inKey);
        switchKeyForTlwe(output, ksKey, input, param);
        auto res = symDecTlweToInt(output, tlweKey, param.torusBase);
        printf("Input plain: %d, Output res: %d.\n", plainMsg, res);
        ASSERT_EQ(plainMsg, res);
    }
    printBanner("KSK");
}