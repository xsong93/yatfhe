//
// Created by Xintong Song on 2024/2/21.
//

#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yautil/numeric_functions.h"
#include "yautil/tool.h"

TEST(PolynomialTest, X) {
    int q = 32;
    int p = 4;
    int N = 32;
    TorusPolynomial v(N);
    std::vector<double> d(N);
    generateTestPolynomial(v, p, 2*N);
//    printArray(v.coeffs, "v");
    for (int i = 0; i < v.N; i++) {
        d[i] = torus32ToDouble(v.coeffs[i]);
    }
    printArray(d, "d");

//    dToT32Test(-3.4);
//    dToT32Test(3.6);
    std::cout << torus32ToDouble(doubleToTorus32(0.51));
    std::cout << torus32ToDouble(doubleToTorus32(-0.49));
//    std::cout << ">>>>>>>> NTT test passed! <<<<<<<<" << std::endl;
}