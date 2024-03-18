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
    YatfheParameters param {};
    DecomposedData decomp(param.ksLevel);
    for (auto t = 0; t < 100; t++) {
        std::vector<Torus> data(5000);
        initCoeffsViaUniformDistribution(data);
        for (auto d : data) {
            signedGadgetDecomposition(decomp, d, param);
            auto recons = recompose(decomp, param);
            ASSERT_EQ(d, recons);
        }
    }
    std::cout << ">>>>>>>> SignedDecomp test passed! <<<<<<<<" << std::endl;
}

TEST(DecomposeTest, DecomposeTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    DecomposedData out(param.ksLevel);
    for (auto i = 0; i < 100; i++) {
        std::vector<Torus> data(5000);
        initCoeffsViaUniformDistribution(data);
        for (auto d : data) {
            gadgetDecompose(out, d, param);
            auto z = recompose(out, param);
            ASSERT_EQ(z, d);
        }
    }
    std::cout << ">>>>>>>> Decompose test passed! <<<<<<<<" << std::endl;
}

TEST(DecomposeOverBTest, DecomposeOverBTest) {
    YatfheParameters param {};
    param.radixBits = 4;
    param.ksLevel = 8;
    auto res = decomposeOverB(1, param);
    printArray(res, "res");
    DecomposedData decomp(param.ksLevel);
    std::vector<Torus> data(5000);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
//        printArray(decomp.value, "decomp");
        int out{0};
        for (auto i = 0; i < res.size(); i++) {
            out += decomp.value[i] * res[i] * decomp.sign;
        }
        ASSERT_EQ(out, d);
//        std::cout << out * decomp.sign << endl;
    }
    std::cout << ">>>>>>>> DecomposeOverB test passed! <<<<<<<<" << std::endl;
}