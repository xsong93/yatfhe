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

TEST(TLWE, ENC_DEC) {
    YatfheParameters param {};
    initYatfhe(param);

    TlweKey tlweKey {param};
    genTlweKey(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    Tlwe input {param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg:" << intModP(plain, param.torusBase) << endl;
    cout << "mu:" << mu << endl;
    cout << "decPre:" << symDecTlweToInt(input, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(plain, param.torusBase), symDecTlweToInt(input, tlweKey, param.torusBase));
    printBanner("TLWE.ENC_DEC");
}

TEST(TLWE, ADD_SUB) {
    YatfheParameters param {};
    param.torusBase = 8;
    initYatfhe(param);

    TlweKey tlweKey {param};
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

    printBanner("TLWE.ADD_SUB");
}

TEST(TLWE, MULT) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = param.torusBits / param.radixBits;
    param.lApprox = param.l - 1;
    initYatfhe(param);

    TlweKey tlweKey {param.n};

    genTlweKey(tlweKey);

    int p1 = 1;
    int p2 = 3;

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

    printBanner("TLWE.MULT");
}