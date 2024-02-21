//
// Created by Xintong Song on 2024/2/21.
//

#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

TEST(PolynomialTest, X) {
    TorusPolynomial v(1024);
    generateTestPolynomial(v, 8, 2048);
    printArray(v.coeffs, "v");
//    std::cout << ">>>>>>>> NTT test passed! <<<<<<<<" << std::endl;
}