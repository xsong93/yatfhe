//
// Created by Xintong Song on 2025/10/28.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

TEST(SERIALIZATION, POLYNOMIAL) {
    YatfheParameters param {};
    param.N = 64;
    initYatfhe(param);

    TorusPolynomial polyT32A{param.N};

    for (size_t i = 0; i < param.N; i++) {
        polyT32A.coeffs[i] = longModP(i, 4);
    }

    printArray(polyT32A.coeffs, "polyT32A");
    std::ostringstream oss(std::ios::binary);
    vector<char> buffer;
    serialize(polyT32A, oss);
    const std::string &data = oss.str();
    buffer.insert(buffer.end(), data.begin(), data.end());
    std::ofstream outFile("SERIALIZATION_TEST_polyT.bin", std::ios::binary | std::ios::trunc);
    outFile.write(buffer.data(), buffer.size());
    buffer.clear();
    outFile.close();

    TorusPolynomial read;
    std::ifstream inFile("SERIALIZATION_TEST_polyT.bin", std::ios::binary);
    deserialize(read, inFile);
    printArray(read.coeffs, "read");

    ASSERT_EQ(read.coeffs, polyT32A.coeffs);

    printBanner("SERIALIZATION.POLYNOMIAL");
}