//
// Created by Xintong Song on 2024/4/29.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

std::array<std::array<int, 4>, 4> KEY_PATTERNS2 = {{
    {1, 0, 0, 0},  // 0b00
    {0, 1, 0, 0},  // 0b01
    {0, 0, 1, 0},  // 0b10
    {0, 0, 0, 1}   // 0b11
}};

std::array<std::array<int, 8>, 8> KEY_PATTERNS3 = {{
    {1, 0, 0, 0, 0, 0, 0, 0},  // 0b000
    {0, 1, 0, 0, 0, 0, 0, 0},  // 0b001
    {0, 0, 1, 0, 0, 0, 0, 0},  // 0b010
    {0, 0, 0, 1, 0, 0, 0, 0},  // 0b011
    {0, 0, 0, 0, 1, 0, 0, 0},  // 0b100
    {0, 0, 0, 0, 0, 1, 0, 0},  // 0b101
    {0, 0, 0, 0, 0, 0, 1, 0},  // 0b110
    {0, 0, 0, 0, 0, 0, 0, 1}   // 0b111
}};

TEST(BOOTSTRAPPING, GROUP2_KEYGEN) {
    YatfheParameters param {};
    param.group = 3;
    param.n = param.group * 7;
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKey bsk {param};
    genBootstrappingKey(bsk, trgswKey, tlweKey, param);

    std::vector<int> idx;
    for (size_t i = 0; i < tlweKey.n; i = i + param.group) {
        int combined = 0;
        for (int i2 = 0; i2 < param.group; i2++) {
            const auto s = tlweKey.s[i + i2] << (param.group - 1 - i2);
            combined |= s;
        }
        idx.push_back(combined);
    }

    std::vector<int> keys;
    for (size_t i = 0; i < bsk.n; i++) {
        keys.push_back(decryptTrgswNtt(bsk.bskDft[i], param, trgswKey));
    }

    printArray(tlweKey.s, "s");
    printArray(idx, "idx");
    printArray(keys, "ks");
    auto batch = 1 << param.group;
    if (param.group == 2) {
        for (size_t i = 0; i < idx.size(); i++) {
            for (size_t j = 0; j < batch; j++) {
                ASSERT_EQ(KEY_PATTERNS2[idx[i]][j], keys[i*batch + j]);
            }
        }
    } else if (param.group == 3) {
        for (size_t i = 0; i < idx.size(); i++) {
            for (size_t j = 0; j < batch; j++) {
                ASSERT_EQ(KEY_PATTERNS3[idx[i]][j], keys[i*batch + j]);
            }
        }
    }
    printBanner("GROUP2_KEYGEN");
}