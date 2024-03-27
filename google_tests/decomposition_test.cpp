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
    std::vector<Torus> data (5000);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
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
        signedGadgetDecomposition(decomp, d, param);
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
    YatfheParameters param{};
    param.radixBits = 4;
    param.ksLevel = 8;
    Trlwe in {2, 4};
    in.a[0].coeffs = {-1,-2,-3,-4};
    in.a[1].coeffs = {2,2,-2,-2};
    in.b.coeffs = {1<<24, 1<<16, 1<<8, 1};
    printTrlweAB(in, "in");
    DecomposedTrlwe out {param.ksLevel, 2, 4};
    gadgetDecomposeTrlwe(out, in, param);
    printDecomposedTrlweAB(out, "out");
}