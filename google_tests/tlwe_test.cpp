//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

TEST(EncDecTest, EncDecTest) {
    const YatfheParameters param {};

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    Tlwe input {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout <<"msg:"<<torus32ToDouble(mu)<<endl;
    cout <<"decPre:"<<symDecTlweSample(input, tlweKey, param.torusBase)<<endl;
    ASSERT_FLOAT_EQ(torus32ToDouble(mu), symDecTlweSample(input, tlweKey, param.torusBase));
    printBanner("EncDecTest");
}

// handles overflow naturally in 32-bits signed int field
TEST(AddSubTest, AddSubTest) {
    YatfheParameters param {};
    param.torusBase = 8;

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    Torus mu1 = doubleToTorus32(3.0 / param.torusBase);
    Torus mu2 = doubleToTorus32(2.0 / param.torusBase);

    Tlwe input1 {param.n};
    Tlwe input2 {param.n};
    Tlwe output {param.n};

    symEncTlweSample(input1, mu1, tlweKey);
    symEncTlweSample(input2, mu2, tlweKey);

    lweAdd(output, input1, input2);
    cout << "plain arithmetic: " << torus32ToDouble(mu1 + mu2) <<endl;
    cout << "dec res: " << symDecTlweSample(output, tlweKey, param.torusBase) << endl;
    ASSERT_FLOAT_EQ(torus32ToDouble(mu1 + mu2), symDecTlweSample(output, tlweKey, param.torusBase));

    lweSub(output, input1, input2);
    cout << "plain arithmetic: " << torus32ToDouble(mu1 - mu2) <<endl;
    cout << "dec res: " << symDecTlweSample(output, tlweKey, param.torusBase) << endl;
    ASSERT_FLOAT_EQ(torus32ToDouble(mu1 - mu2), symDecTlweSample(output, tlweKey, param.torusBase));

    printBanner("AddSubTest");
}