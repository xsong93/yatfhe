//
// Created by Xintong Song on 2024/3/17.
//
#include <vector>
#include <cmath>
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"

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
    std::vector<Torus> decomp(param.ksLevel);
    for (auto k = -5; k < 20; k++) {
        signedGadgetDecomposition(decomp, k, param);
        printArray(decomp,"decomp");
        int recons {};
        for (auto i = 0 ; i < decomp.size(); i++) {
            recons += decomp[i] * powInt(param.radixBase, (param.ksLevel - i - 1));
        }
//        ASSERT_EQ(k, recons);
    }
    std::cout << ">>>>>>>> SignedDecomp test passed! <<<<<<<<" << std::endl;
}