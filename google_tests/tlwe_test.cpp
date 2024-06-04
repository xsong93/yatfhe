//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

TEST(TlweTest, EncDecTest) {
    const YatfheParameters param {};

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    int plain = genIntUniformDist(-4, 3);
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    Tlwe input {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout << "msg:" << intModP(plain, param.torusBase) << endl;
    cout << "decPre:" << symDecTlweSampleToInt(input, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(plain, param.torusBase), symDecTlweSampleToInt(input, tlweKey, param.torusBase));
    printBanner("EncDecTest");
}
//
//TEST(TlweTest, EncDecTest2) {
//    const YatfheParameters param {};
//
//    TlweKey tlweKey {param.n, 0};
//    lweKeyGen(tlweKey);
//
//    int mu = 0;
//    Tlwe input {param.n};
//    symEncTlweSample(input, mu, tlweKey);
//    printTlweAB(input, "in32");
//
//    Tlwe input8 {param.n};
//    for (auto i = 0; i < param.n; i++) {
//        input8.a[i] = modSwitchFromTorus32(input.a[i], 128);
//    }
//    input8.b = modSwitchFromTorus32(input.b, 128);
//    printTlweAB(input8, "in8");
//
//    Torus aXs = 0;
//    for (auto i = 0; i < tlweKey.n; i++) {
//        aXs += input8.a[i] * tlweKey.s[i];
//        aXs = aXs % (128);
//    }
//    cout <<"aXs:"<<aXs<<endl;
//    cout <<"b-aXs:"<<(input8.b - aXs) % (128)<<endl;
//
//    cout <<"msg:"<<mu<<endl;
//    cout <<"decPre:"<<symDecTlweSampleToInt(input, tlweKey, 128)<<endl;
//    ASSERT_EQ(mu, symDecTlweSampleToInt(input, tlweKey, 128));
//    printBanner("EncDecTest2");
//}

// handles overflow naturally in 32-bits signed int field
TEST(TlweTest, AddSubTest) {
    YatfheParameters param {};
    param.torusBase = 8;

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    Torus mu1 = modSwitchToTorus32(3, param.torusBase);
    Torus mu2 = modSwitchToTorus32(2, param.torusBase);

    Tlwe input1 {param.n};
    Tlwe input2 {param.n};
    Tlwe output {param.n};

    symEncTlweSample(input1, mu1, tlweKey);
    symEncTlweSample(input2, mu2, tlweKey);

    lweAdd(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorus32(mu1 + mu2, param.torusBase) <<endl;
    cout << "dec res: " << symDecTlweSampleToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorus32(mu1 + mu2, param.torusBase), symDecTlweSampleToInt(output, tlweKey, param.torusBase));

    lweSub(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorus32(mu1 - mu2, param.torusBase) <<endl;
    cout << "dec res: " << symDecTlweSampleToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorus32(mu1 - mu2, param.torusBase), symDecTlweSampleToInt(output, tlweKey, param.torusBase));

    printBanner("AddSubTest");
}