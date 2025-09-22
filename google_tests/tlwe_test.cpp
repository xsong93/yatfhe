//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yatfhe/tlev.h"
#include "yautil/initializer.h"

TEST(TlweTest, RED_TEST) {
    YatfheParameters param {};
    initYatfhe(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    Tlwe input {param.n};
    symEncTlwe(input, mu, tlweKey);

    // saw off
    int thres = 16;
    for (size_t j = 0; j < param.N; j++) {
        auto tmp = input.a[j] >> thres;
        input.a[j] = tmp << thres;
    }
    auto tmp = input.b >> thres;
    input.b = tmp << thres;

    cout << "msg:" << intModP(plain, param.torusBase) << endl;
    cout << "mu:" << mu << endl;
    cout << "decPre:" << symDecTlweToInt(input, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(plain, param.torusBase), symDecTlweToInt(input, tlweKey, param.torusBase));
    printBanner("RED_TEST");
}

TEST(TlweTest, EncDecTest) {
    YatfheParameters param {};
    initYatfhe(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
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

    Torus mu1 = modSwitchToTorus32(3, param.torusBase);
    Torus mu2 = modSwitchToTorus32(2, param.torusBase);

    Tlwe input1 {param.n};
    Tlwe input2 {param.n};
    Tlwe output {param.n};

    symEncTlwe(input1, mu1, tlweKey);
    symEncTlwe(input2, mu2, tlweKey);

    addTlwe(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorus32(mu1 + mu2, param.torusBase) <<endl;
    cout << "dec res: " << symDecTlweToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorus32(mu1 + mu2, param.torusBase), symDecTlweToInt(output, tlweKey, param.torusBase));

    subTlwe(output, input1, input2);
    cout << "plain arithmetic: " << modSwitchFromTorus32(mu1 - mu2, param.torusBase) <<endl;
    cout << "dec res: " << symDecTlweToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(modSwitchFromTorus32(mu1 - mu2, param.torusBase), symDecTlweToInt(output, tlweKey, param.torusBase));

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

    Torus mu1 = modSwitchToTorus32(p1, param.torusBase);

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

//todo
TEST(TlweTest, todotest) {
    YatfheParameters param {};
    param.torusBase = 8;
    initYatfhe(param);
    int scale = Q_32 >> 1;

    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);

    vector<Tlwe> tlwes(param.n, Tlwe{param.n});
    for (size_t i = 0; i < param.n; i++) {
        symEncTlwe(tlwes[i], modSwitchToTorus32(tlweKey.s[i], scale), tlweKey);
        cout << "s:" << symDecTlweToInt(tlwes[i], tlweKey, scale) << ", ";
    }
    cout << endl;
    Tlwe sample{param.n};
    int pt = 2;
    auto mu = modSwitchToTorus32(pt, param.torusBase);
    symEncTlwe(sample, mu, tlweKey);

    auto dec = symDecTlweToInt(sample, tlweKey, param.torusBase);
    cout << "dec:" << dec << endl;

    ScaledTlwe inputModN2 {scale, param.n};
    rescaleTlweFromTorus32(inputModN2, sample);

    vector<Integer> redA(param.n);
    for (size_t i = 0; i < param.n; i++) {
        redA[i] = inputModN2.a[i] % scale;
    }
    Integer redB = inputModN2.b % scale;
    printArray(redA, "redA");
    cout << "redB:" << redB << endl;

    Tlwe zero{param.n};
    symEncTlwe(zero, 0, tlweKey);
    zero.b = addTorus(zero.b, modSwitchToTorus32(redB, scale));
    cout << "decB:" << symDecTlweToInt(zero, tlweKey, scale) << endl;

    for (size_t i = 0; i < param.n; i++) {
        for (size_t j = 0; j < param.n; j++) {
            tlwes[i].a[j] = multTorus(tlwes[i].a[j], redA[i]);
        }
        tlwes[i].b = multTorus(tlwes[i].b, redA[i]);
        cout << "as:" << symDecTlweToInt(tlwes[i], tlweKey, scale) << ", ";
        for (size_t j = 0; j < param.n; j++) {
            zero.a[j] = subTorus(zero.a[j], tlwes[i].a[j]);
        }
        zero.b = subTorus(zero.b, tlwes[i].b);
    }
    cout << endl;

    auto dec1 = symDecTlweToInt(zero, tlweKey, 8);
    cout << "dec1:" << dec1 << endl;
    cout <<"err0:" << calTlweError(sample, tlweKey, mu) << endl;
    cout <<"err1:" << calTlweError(zero, tlweKey, modSwitchToTorus32(pt, 8)) << endl;

    printBanner("TlweMultTest");
}