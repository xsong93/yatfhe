//
// Created by Xintong Song on 2024/3/20.
//
#include <vector>
#include "gtest/gtest.h"
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

TEST(PolyRotTest, PolyRotTest) {
    int N = 4;
    std::vector<Torus> m7 {-4,1,2,3};
    std::vector<Torus> m6 {-3,-4,1,2};
    std::vector<Torus> m5 {-2,-3,-4,1};
    std::vector<Torus> m4 {-1,-2,-3,-4};
    std::vector<Torus> m3 {4,-1,-2,-3};
    std::vector<Torus> m2 {3,4,-1,-2};
    std::vector<Torus> m1 {2,3,4,-1};
    std::vector<Torus> zo {1,2,3,4};
    std::vector<std::vector<Torus>> dic (N * 2, std::vector<Torus> (N));
    dic[0] = zo; dic[8] = zo; dic[16] = zo;
    dic[1] = m7; dic[9] = m7;
    dic[2] = m6; dic[10] = m6;
    dic[3] = m5; dic[11] = m5;
    dic[4] = m4; dic[12] = m4;
    dic[5] = m3; dic[13] = m3;
    dic[6] = m2; dic[14] = m2;
    dic[7] = m1; dic[15] = m1;
    TorusPolynomial a {N};
    TorusPolynomial b {N};
    for (auto i = 0; i < N ; i++) {
        a.coeffs[i] = i + 1;
    }
    printArray(a.coeffs, "a");
    for (auto i = - N * 2; i <= N * 2; i++) {
        torusPolynomialRotate(b, i, a);
        printArray(b.coeffs, "b " + to_string(i));
        ASSERT_EQ(b.coeffs, dic[i + N * 2]);
    }
}