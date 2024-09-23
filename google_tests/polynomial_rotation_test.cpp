//
// Created by Xintong Song on 2024/3/20.
//
#include <vector>
#include "gtest/gtest.h"
#include "yatfhe/torus.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

std::vector<Torus> vecSub(const std::vector<Torus>& a, const std::vector<Torus>& b) {
    std::vector<Torus> tmp (a.size());
    for (auto i = 0; i < a.size(); i++) {
        tmp[i] = a[i] - b[i];
    }
    return tmp;
}

std::vector<int8_t> vecSub8(const std::vector<int8_t>& a, const std::vector<int8_t>& b) {
    std::vector<int8_t> tmp (a.size());
    for (size_t i = 0; i < a.size(); i++) {
        tmp[i] = a[i] - b[i];
    }
    return tmp;
}

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
    std::vector<std::vector<Torus>> rotatedVec (N * 4 + 1, std::vector<Torus> (N));
    rotatedVec[0] = zo; rotatedVec[8] = zo; rotatedVec[16] = zo;
    rotatedVec[1] = m7; rotatedVec[9] = m7;
    rotatedVec[2] = m6; rotatedVec[10] = m6;
    rotatedVec[3] = m5; rotatedVec[11] = m5;
    rotatedVec[4] = m4; rotatedVec[12] = m4;
    rotatedVec[5] = m3; rotatedVec[13] = m3;
    rotatedVec[6] = m2; rotatedVec[14] = m2;
    rotatedVec[7] = m1; rotatedVec[15] = m1;
    TorusPolynomial a {N};
    TorusPolynomial b {N};
    for (auto i = 0; i < N ; i++) {
        a.coeffs[i] = i + 1;
    }
    printArray(a.coeffs, "a");
    for (auto i = - N * 2; i <= N * 2; i++) {
        torusPolynomialRotate(b, i, a);
        printArray(b.coeffs, "b " + to_string(i));
        ASSERT_EQ(b.coeffs, rotatedVec[i + N * 2]);
    }
    for (auto i = - N * 2; i <= N * 2; i++) {
        torusPolynomialRotateMinusOne(b, i, a);
        printArray(b.coeffs, "b-1 " + to_string(i));
        ASSERT_EQ(b.coeffs, vecSub(rotatedVec[i + N * 2], a.coeffs));
    }
    printBanner("PolyRotTest");
}

TEST(Poly8RotTest, Poly8RotTest) {
    int N = 4;
    std::vector<int8_t> m7 {-4,1,2,3};
    std::vector<int8_t> m6 {-3,-4,1,2};
    std::vector<int8_t> m5 {-2,-3,-4,1};
    std::vector<int8_t> m4 {-1,-2,-3,-4};
    std::vector<int8_t> m3 {4,-1,-2,-3};
    std::vector<int8_t> m2 {3,4,-1,-2};
    std::vector<int8_t> m1 {2,3,4,-1};
    std::vector<int8_t> zo {1,2,3,4};
    std::vector<std::vector<int8_t>> rotatedVec (N * 4 + 1, std::vector<int8_t> (N));
    rotatedVec[0] = zo; rotatedVec[8] = zo; rotatedVec[16] = zo;
    rotatedVec[1] = m7; rotatedVec[9] = m7;
    rotatedVec[2] = m6; rotatedVec[10] = m6;
    rotatedVec[3] = m5; rotatedVec[11] = m5;
    rotatedVec[4] = m4; rotatedVec[12] = m4;
    rotatedVec[5] = m3; rotatedVec[13] = m3;
    rotatedVec[6] = m2; rotatedVec[14] = m2;
    rotatedVec[7] = m1; rotatedVec[15] = m1;
    Int8Polynomial a {N};
    Int8Polynomial b {N};
    for (auto i = 0; i < N ; i++) {
        a.coeffs[i] = i + 1;
    }
    printArray(a.coeffs, "a");
    for (auto i = - N * 2; i <= N * 2; i++) {
        int8PolynomialRotate(b, i, a, 256);
        printArray(b.coeffs, "b " + to_string(i));
        ASSERT_EQ(b.coeffs, rotatedVec[i + N * 2]);
    }
    for (auto i = - N * 2; i <= N * 2; i++) {
        int8PolynomialRotateMinusOne(b, i, a, 256);
        printArray(b.coeffs, "b-1 " + to_string(i));
        ASSERT_EQ(b.coeffs, vecSub8(rotatedVec[i + N * 2], a.coeffs));
    }

    for (auto i = 0; i < N ; i++) {
        a.coeffs[i] = -128 - i;
    }
    for (auto i = - N * 2; i <= N * 2; i++) {
        int8PolynomialRotate(b, i, a, 256);
        printArray(b.coeffs, "-128-1 " + to_string(i));
//        ASSERT_EQ(b.coeffs, vecSub8(rotatedVec[i + N * 2], a.coeffs));
    }
    printBanner("Poly8RotTest");
}