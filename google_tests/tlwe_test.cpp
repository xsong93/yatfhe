//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"

TEST(EncDecTest, EncDecTest) {
    const YatfheParameters param {};

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey, param.n);

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    Tlwe input {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout <<"msg:"<<torus32ToDouble(mu)<<endl;
    cout <<"decPre:"<<symDecTlweSample(input, tlweKey, param.torusBase)<<endl;
    ASSERT_FLOAT_EQ(torus32ToDouble(mu), symDecTlweSample(input, tlweKey, param.torusBase));
}