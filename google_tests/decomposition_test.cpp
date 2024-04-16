//
// Created by Xintong Song on 2024/3/17.
//
#include <vector>
#include <cmath>
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"

UnsignedInteger powInt(UnsignedInteger base, UnsignedInteger exponent) {
    UnsignedInteger result = 1;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    return result;
}

TEST(SignedDecompTest, SignedDecompTest) {
    YatfheParameters param {};
    DecomposedData decomp {param.ksLevel};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data);
//    data = {-10,-9,9,10};
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
        printf("in: %d, ", d);
        printArray(decomp.value, "decomp");
        auto recons = selfRecompose(decomp, param);
        ASSERT_EQ(d, recons);
    }
    printBanner("SignedDecomp");
}

TEST(DecomposeTest, DecomposeTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    DecomposedData out {param.ksLevel};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
        printf("in: %d, ", d);
        printArray(out.value, "decomp");
        auto z = selfRecompose(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("Decompose");
}

TEST(DecomposeOverBSingleStage, DecomposeOverBSingleStage) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    int mult = 5;
    std::vector<Integer> rhs(param.ksLevel);
    decomposeOverB(rhs, mult, param);
    printArray(rhs, to_string(mult) + " decomposeOverB");
    DecomposedData decomp {param.ksLevel};
    std::vector<Torus> data(10);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
//        signedGadgetDecomposition(decomp, d, param); // both correct
        gadgetDecompose(decomp, d, param); // both correct
        printf("in: %d, ", d);
        printArray(decomp.value, "decomp");
        auto out = recomposeTwoParts(decomp, rhs);
        ASSERT_EQ(out, d * mult);
    }
    printBanner("DecomposeOverBSingleStage");
}

TEST(DecomposeOverBMultiStages, DecomposeOverBMultiStages) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    int mult = -5;

    // first decomp
    std::vector<Integer> rhs(param.ksLevel);
    decomposeOverB(rhs, mult, param);
    printArray(rhs, to_string(mult) + " decomposeOverB");

    // second decomp
    vector<DecomposedData> decompRhs {rhs.size(), DecomposedData(param.ksLevel)};
    for (auto i = 0; i < rhs.size(); i++) {
        gadgetDecompose(decompRhs[i], rhs[i], param);
        printArray(decompRhs[i].value, "rhs level " + to_string(i));
    }

    // recomp second decomposed data test
    vector<Integer> recompRhs(rhs.size());
    for (auto i = 0; i < decompRhs.size(); i++) {
        recompRhs[i] = selfRecompose(decompRhs[i], param);
        ASSERT_EQ(recompRhs[i], rhs[i]);
    }
    printArray(recompRhs, "recompRhs");

    DecomposedData decompL1 {param.ksLevel};
    DecomposedData recompL1 {param.ksLevel};
    Torus data = genIntUniformDist(TorusMin, TorusMax);

//  signedGadgetDecomposition(decomp, data, param); // both correct
    gadgetDecompose(decompL1, data, param); // both correct
    printf("in: %d, ", data);
    printArray(decompL1.value, "decomp");

    // recomp first level
    recomposeFirstHalf(recompL1, decompL1, decompRhs);
    printArray(recompL1.value, "recompL1");

    // recomp second level
    auto out = selfRecompose(recompL1, param); // equivalent to recomposeTwoParts(recompL1, decompOneOverR)
    printf("out = %d, data * mult = %d\n", out, data * mult);
    ASSERT_EQ(out, data * mult);

    printBanner("DecomposeOverBMultiStages");
}

TEST(DecomposeTrlweTest, DecomposeTrlweTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.l = 8;
    param.k = 2;
    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};
    TrlweDft recompDft {param.k, param.N};

//    param.N = 4;
//    in.a[0].coeffs = {-1,-2,-3,-4};
//    in.a[1].coeffs = {2,2,-2,-2};
//    in.b.coeffs = {1<<24, 1<<16, 1<<8, 1};
//    applyNttForAB(inDft, in);
//    printTrlweAB(in, "in");

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    TrlweKey trlweKey {param.k, param.N};
    trlweKeyGen(trlweKey);
    symEncTrlweSingleSample(in, inDft, trlweKey, mu, param.lweStdDev);
    Trlwe intt {param.k, param.N};
    applyInttForAB(intt, inDft);

    // decompose
    DecomposedTrlwe out {param};
    gadgetDecomposeTrlwe(out, in, param);
    gadgetDecomposeTrlweNtt(out, inDft, param);
//    printDecomposedTrlweAB(out, "out");
//    printDecomposedTrlweNttAB(out, "outNtt");

    // recompose original
    recomposeTrlwe(recomp, out, param);
    printTrlweAB(in, "original");
    printTrlweAB(intt, "intt");
    printTrlweAB(recomp, "recomp");
    for (auto j = 0; j < in.b.N; j++) {
        for (auto i = 0 ; i < in.k; i++) {
            ASSERT_EQ(in.a[i].coeffs[j], intt.a[i].coeffs[j]);
            ASSERT_EQ(in.a[i].coeffs[j], recomp.a[i].coeffs[j]);
        }
        ASSERT_EQ(in.b.coeffs[j], intt.b.coeffs[j]);
        ASSERT_EQ(in.b.coeffs[j], recomp.b.coeffs[j]);
    }

    // todo: test decompose/ntt order
    // recompose ntt
    recomposeTrlweNtt(recompDft, out, param);
    printTrlweDftAB(inDft, "inDft");
    printTrlweDftAB(recompDft, "recompDft");
    applyInttForAB(intt, recompDft);
    printTrlweAB(intt, "recompDft:intt");

    for (auto j = 0; j < in.b.N; j++) {
        for (auto i = 0 ; i < in.k; i++) {
            ASSERT_EQ(in.a[i].coeffs[j], intt.a[i].coeffs[j]);
        }
        ASSERT_EQ(in.b.coeffs[j], intt.b.coeffs[j]);
    }
    printBanner("DecomposeTrlweTest");
}

TEST(DecomposedAddSub, DecomposedAddSub) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    DecomposedData da {param.ksLevel};
    DecomposedData db {param.ksLevel};
    DecomposedData dr {param.ksLevel};
    int a;
    int b;
    int t = 5000;
    while (t-- > 0) {
        a = genIntUniformDist(INT_MIN, INT_MAX);
        b = genIntUniformDist(INT_MIN, INT_MAX);
        gadgetDecompose(da, a, param);
        gadgetDecompose(db, b, param);
        for (auto i = 0; i < param.ksLevel; i++) {
            dr.value[i] = da.value[i] * da.sign + db.value[i] * db.sign;
        }
        auto z = selfRecompose(dr, param);
        printf("a + b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a + b);
        for (auto i = 0; i < param.ksLevel; i++) {
            dr.value[i] = da.value[i] * da.sign - db.value[i] * db.sign;
        }
        z = selfRecompose(dr, param);
        printf("a - b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a - b);
    }
    printBanner("DecomposedAddSub");
}