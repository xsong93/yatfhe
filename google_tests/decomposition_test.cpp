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
        auto recons = recompose(decomp, param);
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
        auto z = recompose(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("Decompose");
}

TEST(DecomposeOverBTest, DecomposeOverBTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    int mult = 5;
    auto rhs = decomposeOverB(mult, param);
    printArray(rhs, "1 decomposeOverB");
    DecomposedData decomp {param.ksLevel};
    std::vector<Torus> data(10);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param); // both correct
//        gadgetDecompose(decomp, d, param); // both correct
        printf("in: %d, ", d);
        printArray(decomp.value, "decomp");
        int out {0};
        for (auto i = 0; i < rhs.size(); i++) {
            out += decomp.value[i] * rhs[i] * decomp.sign;
        }
        ASSERT_EQ(out, d * mult);
    }
    printBanner("DecomposeOverB");
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
    DecomposedTrlwe out {param.l, param.k, param.N};
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