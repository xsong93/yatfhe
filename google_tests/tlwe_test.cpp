//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yatfhe/tglev.h"
#include "yautil/initializer.h"

TEST(TlweTest, RED_TEST) {
    YatfheParameters param {};
    yatfheInit(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    Tlwe input {param.n};
    symEncTlweSample(input, mu, tlweKey);

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
    cout << "decPre:" << symDecTlweSampleToInt(input, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(plain, param.torusBase), symDecTlweSampleToInt(input, tlweKey, param.torusBase));
    printBanner("RED_TEST");
}

TEST(TlweTest, EncDecTest) {
    YatfheParameters param {};
    yatfheInit(param);

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    int plain = genIntUniformDist(-param.torusBase/2, param.torusBase/2 - 1);
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    Tlwe input {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout << "msg:" << intModP(plain, param.torusBase) << endl;
    cout << "mu:" << mu << endl;
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
    yatfheInit(param);

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

TEST(TlweTest, MultTest) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.l = 4;
    yatfheInit(param);

    TlweKey tlweKey {param.n, param.lweStdDev};

    lweKeyGen(tlweKey);

    int p1 = 1;
    int p2 = 2047;

    Torus mu1 = modSwitchToTorus32(p1, param.torusBase);

    Tglev tglev {param.l, param.n};
    Tlwe output {param.n};

    tglevEnc(tglev, tlweKey, mu1, param);

//    for (auto l = 0; l < tglev.l; l++) {
//        printArray(tglev.tlwes[l].a, "a");
//        cout << "b: " << tglev.tlwes[l].b << endl;
//    }

    tglevMultConst(output, tglev, p2, param);
    cout << "p1 * p2: " << intModP(p1 * p2, param.torusBase) <<endl;
    cout << "c1 * p2: " << symDecTlweSampleToInt(output, tlweKey, param.torusBase) << endl;
    ASSERT_EQ(intModP(p1 * p2, param.torusBase), symDecTlweSampleToInt(output, tlweKey, param.torusBase));

    printBanner("TlweMultTest");
}

//todo
TEST(TlweTest, todotest) {
    YatfheParameters param {};
    param.torusBase = 8;
    yatfheInit(param);
    int scale = Q_32 >> 1;

    TlweKey tlweKey {param.n, param.lweStdDev};
    lweKeyGen(tlweKey);

    vector<Tlwe> tlwes(param.n, Tlwe{param.n});
    for (size_t i = 0; i < param.n; i++) {
        symEncTlweSample(tlwes[i], modSwitchToTorus32(tlweKey.s[i], scale), tlweKey);
        cout << "s:" << symDecTlweSampleToInt(tlwes[i], tlweKey, scale) << ", ";
    }
    cout << endl;
    Tlwe sample{param.n};
    int pt = 2;
    auto mu = modSwitchToTorus32(pt, param.torusBase);
    symEncTlweSample(sample, mu, tlweKey);

    auto dec = symDecTlweSampleToInt(sample, tlweKey, param.torusBase);
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
    symEncTlweSample(zero, 0, tlweKey);
    zero.b = modAddT32(zero.b, modSwitchToTorus32(redB, scale));
    cout << "decB:" << symDecTlweSampleToInt(zero, tlweKey, scale) << endl;

    for (size_t i = 0; i < param.n; i++) {
        for (size_t j = 0; j < param.n; j++) {
            tlwes[i].a[j] = modMulT32(tlwes[i].a[j], redA[i]);
        }
        tlwes[i].b = modMulT32(tlwes[i].b, redA[i]);
        cout << "as:" << symDecTlweSampleToInt(tlwes[i], tlweKey, scale) << ", ";
        for (size_t j = 0; j < param.n; j++) {
            zero.a[j] = modSubT32(zero.a[j], tlwes[i].a[j]);
        }
        zero.b = modSubT32(zero.b, tlwes[i].b);
    }
    cout << endl;

    auto dec1 = symDecTlweSampleToInt(zero, tlweKey, 8);
    cout << "dec1:" << dec1 << endl;
    cout <<"err0:" << calTlweError(sample, tlweKey, mu) << endl;
    cout <<"err1:" << calTlweError(zero, tlweKey, modSwitchToTorus32(pt, 8)) << endl;

    printBanner("TlweMultTest");
}