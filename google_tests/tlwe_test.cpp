//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yatfhe/tlev.h"
#include "yautil/initializer.h"

TEST(TlweTest, EncDecTest) {
    YatfheParameters param {};
    initYatfhe(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    Tlwe input {param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg:" << intModP(plain, param.torusBase) << endl;
    cout << "mu:" << mu << endl;
    cout << "decPre:" << symDecTlweToInt(input, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(plain, param.torusBase), symDecTlweToInt(input, tlweKey, param.torusBase));
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
//        input8.a[i] = modSwitchFromTorus32(input.a[i], 1024);
//    }
//    input8.b = modSwitchFromTorus32(input.b, 1024);
//    printTlweAB(input8, "in8");
//
//    Torus aXs = 0;
//    for (auto i = 0; i < tlweKey.n; i++) {
//        aXs += input8.a[i] * tlweKey.s[i];
//        aXs = aXs % (1024);
//    }
//    cout <<"aXs:"<<aXs<<endl;
//    cout <<"b-aXs:"<<(input8.b - aXs) % (1024)<<endl;
//
//    cout <<"msg:"<<mu<<endl;
//    cout <<"decPre:"<<symDecTlweSampleToInt(input, tlweKey, 1024)<<endl;
//    ASSERT_EQ(mu, symDecTlweSampleToInt(input, tlweKey, 1024));
//    printBanner("EncDecTest2");
//}

// handles overflow naturally in 32-bits signed int field
TEST(TlweTest, AddSubTest) {
    YatfheParameters param {};
    param.torusBase = 8;
    initYatfhe(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);

    Torus mu1 = modSwitchToTorusGeneral(3, param.torusBase, LWE_Q);
    Torus mu2 = modSwitchToTorusGeneral(2, param.torusBase, LWE_Q);

    Tlwe input1 {param.n};
    Tlwe input2 {param.n};
    Tlwe output {param.n};

    symEncTlwe(input1, mu1, tlweKey);
    symEncTlwe(input2, mu2, tlweKey);

    addTlwe(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorusGeneral(addTorus(LWE_Q, mu1, mu2), param.torusBase, LWE_Q) <<endl;
    cout << "dec res: " << symDecTlweToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorusGeneral(addTorus(LWE_Q, mu1, mu2), param.torusBase, LWE_Q), symDecTlweToInt(output, tlweKey, param.torusBase));

    subTlwe(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorusGeneral(subTorus(LWE_Q, mu1, mu2), param.torusBase, LWE_Q) <<endl;
    cout << "dec res: " << symDecTlweToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorusGeneral(subTorus(LWE_Q, mu1, mu2), param.torusBase, LWE_Q), symDecTlweToInt(output, tlweKey, param.torusBase));

    printBanner("AddSubTest");
}

TEST(TlweTest, MultTest) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.l = 4;
    initYatfhe(param);

    TlweKey tlweKey {param.n, param.lweStdDev};

    genTlweKey(tlweKey);

    int p1 = 1;
    int p2 = 2047;

    Torus mu1 = modSwitchToTorusGeneral(p1, param.torusBase, LWE_Q);

    Tlev tglev {param.l, param.n};
    Tlwe output {param.n};

    encTglev(tglev, tlweKey, mu1, param);

//    for (auto l = 0; l < tglev.l; l++) {
//        printArray(tglev.tlwes[l].a, "a");
//        cout << "b: " << tglev.tlwes[l].b << endl;
//    }

    multTglevWithConst(output, tglev, p2, param);
    cout << "p1 * p2: " << intModP(p1 * p2, param.torusBase) <<endl;
    cout << "c1 * p2: " << symDecTlweToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(p1 * p2, param.torusBase), symDecTlweToInt(output, tlweKey, param.torusBase));

    printBanner("TlweMultTest");
}