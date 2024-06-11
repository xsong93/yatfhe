//
// Created by Xintong Song on 2024/6/11.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

TEST(ModSwitchTest, ModDownTest) {
    const YatfheParameters param {};
    const int modPQ = 256;

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    // data gen
    int mu = 1;
    Torus mt = modSwitchToTorus32(mu, 8);
    Tlwe input {param.n};
    symEncTlweSample(input, mt, tlweKey);
    int mu2 = 2;
    Torus mt2 = modSwitchToTorus32(mu2, 8);
    Tlwe input2 {param.n};
    symEncTlweSample(input2, mt2, tlweKey);

    // mod down
    for (auto i = 0; i < param.n; i++) {
        input.a[i] = modSwitchFromTorus32(input.a[i], modPQ);
    }
    input.b = modSwitchFromTorus32(input.b, modPQ);
    for (auto i = 0; i < param.n; i++) {
        input2.a[i] = modSwitchFromTorus32(input2.a[i], modPQ);
    }
    input2.b = modSwitchFromTorus32(input2.b, modPQ);

    // mod add
    Tlwe add {param.n};
    for (auto i = 0; i < add.n; i++) {
        add.a[i] = intModP(input.a[i] + input2.a[i], modPQ);
    }
    add.b = intModP(input.b + input2.b, modPQ);

    // mod up
    for (auto i = 0; i < param.n; i++) {
        add.a[i] = modSwitchToTorus32(add.a[i], modPQ);
    }
    add.b = modSwitchToTorus32(add.b, modPQ);

    cout <<"msg:"<<intModP(mu + mu2, 8)<<endl;
    cout <<"dec:"<<symDecTlweSampleToInt(add, tlweKey, 8)<<endl;
    ASSERT_EQ(intModP(mu + mu2, 8), symDecTlweSampleToInt(add, tlweKey, 8));
    printBanner("ModDownTest");
}