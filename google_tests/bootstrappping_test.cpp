//
// Created by Xintong Song on 2024/4/29.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/key_patterns.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(BOOTSTRAPPING, GROUP2_KEYGEN) {
    YatfheParameters param {};
    param.n = 8;
    param.group = 2;
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
    for (size_t i = 0; i < tlweKey.n; i = i + 2) {
        auto s1 = tlweKey.s[i];
        auto s2 = tlweKey.s[i + 1];
        idx.push_back((s1 << 1) | s2);
    }

    std::vector<int> keys;
    for (size_t i = 0; i < bsk.n; i++) {
        keys.push_back(decryptTrgswNtt(bsk.bskDft[i], param, trgswKey));
    }

    auto batch = 1 << param.group;
    for (size_t i = 0; i < idx.size(); i++) {
        for (size_t j = 0; j < batch; j++) {
            ASSERT_EQ(KEY_PATTERNS2[idx[i]][j], keys[i*batch + j]);
        }
    }
    printArray(tlweKey.s, "s");
    printArray(keys, "ks");
    printBanner("GROUP_KEYGEN");
}