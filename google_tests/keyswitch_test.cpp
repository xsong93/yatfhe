//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(KSK, KSK) {
    YatfheParameters param {};
    param.torusBase = 1 << 3;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    TlweKey tlweKey {param};
    TrlweKey trlweKey {param};
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TlweKey inKey(param.k * param.N, param.rlweNoiseB);
    convertTrlweKeyToTlweKey(inKey, trlweKey);
    for (int i = - param.torusBase / 2; i < param.torusBase / 2; i++) {
        int plainMsg = i;
        Torus mu = modSwitchToTorusGeneral(plainMsg, param.torusBase, LWE_Q);
        Tlwe input {param.N * param.k};
        Tlwe output {param.n};
        symEncTlwe(input, mu, inKey);
        switchKeyForTlwe(output, ksKey, input, param);
        auto res = symDecTlweToInt(output, tlweKey, param.torusBase);
        printf("Input plain: %d, Output res: %d.\n", plainMsg, res);
        ASSERT_EQ(plainMsg, res);
    }
    printBanner("KSK.KSK");
}

//todo
TEST(KSK, PRIVATE_KSK) {
    YatfheParameters param{};
    param.torusBase = 1 << 3;
    param.driftPhase = param.N / param.torusBase / 2;
    initYatfhe(param);

    TlweKey tlweKey{param};
    TrlweKey trlweKey {param};
    PrivateKeySwitchingKey psk{param, param.l};

    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    genPrivateKeySwitchingKey(psk, trlweKey, tlweKey, param);

    for (int i = - param.torusBase / 2; i < param.torusBase / 2; i++) {
        Integer plainMsg = i;
        Torus mu = modSwitchToTorusGeneral(plainMsg, param.torusBase, LWE_Q);
        Tlwe input{param.n}, tmp {param.n};
        Trlwe output{param};
        TorusPolynomial tp{param.N};
        symEncTlwe(input, mu, tlweKey);
        tlweToTrlwePrivateKeySwitching(output, input, psk, param);
        symDecTrlweToInt(tp, output, trlweKey, param.torusBase);
//        printArray(tp.coeffs, "tp");
        printf("Input plain: %d, Output res: %d.\n", plainMsg, tp.coeffs[0]);
        ASSERT_EQ(plainMsg, tp.coeffs[0]);
    }

    printBanner("KSK.PRIVATE_KSK");
}
