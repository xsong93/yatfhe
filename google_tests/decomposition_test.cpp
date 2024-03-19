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

u_int32_t powInt(u_int32_t base, u_int32_t exponent) {
    u_int32_t result = 1;
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
    YatfheParameters param{};
    DecomposedData decomp{param.ksLevel};
    std::vector<Torus> data(5000);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
        auto recons = recompose(decomp, param);
        ASSERT_EQ(d, recons);
    }
    printBanner("SignedDecomp");
}

TEST(DecomposeTest, DecomposeTest) {
    YatfheParameters param{};
    param.radixBits = 4;
    param.ksLevel = 8;
    DecomposedData out{param.ksLevel};
    std::vector<Torus> data(5000);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
        auto z = recompose(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("Decompose");
}

TEST(DecomposeOverBTest, DecomposeOverBTest) {
    YatfheParameters param{};
    param.radixBits = 4;
    param.ksLevel = 8;
    auto rhs = decomposeOverB(1, param);
    printArray(rhs, "1 decomposeOverB");
    DecomposedData decomp{param.ksLevel};
    std::vector<Torus> data(5000);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
//        printArray(decomp.value, "decomp");
        int out{0};
        for (auto i = 0; i < rhs.size(); i++) {
            out += decomp.value[i] * rhs[i] * decomp.sign;
        }
        ASSERT_EQ(out, d);
//        std::cout << out * decomp.sign << endl;
    }
    printBanner("DecomposeOverB");
}