//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

TEST(KSKTest, KSKTest) {
    YatfheParameters param{};
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param.N * param.k, param.n, param.ksLevel};
    lweKeyGen(tlweKey, param.n);
    trlweKeyGen(trlweKey, param.N, param.k);
    tlweKeySwitchingKeyGen(ksKey, trlweKey, tlweKey, param);
    TlweKey inKey(param.k * param.N);
    convertTrlweKeyToTlweKey(inKey, trlweKey);
    param.torusBase = 1 << 8;
    for (int i = - param.torusBase / 2 + 1; i < param.torusBase / 2; i++) {
        double plainMsg = (double)i / param.torusBase;
        Torus mu = doubleToTorus32(plainMsg);
        Tlwe input {param.N * param.k};
        Tlwe output {param.n};
        symEncTlweSample(input, mu, inKey);
        tlweKeySwitch(output, ksKey, input, param);
        auto res = symDecTlweSample(output, tlweKey, param.torusBase);
        printf("Input plain: %f, Output res: %f.\n", plainMsg, res);
        ASSERT_EQ(plainMsg * param.torusBase, res * param.torusBase);
    }
    printBanner("KSK");
}