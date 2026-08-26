//
// Created by Xintong Song on 2025/10/28.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yautil/ya_serializer.h"

void verifyTrlwe(const Trlwe& a, const Trlwe& b) {
    ASSERT_EQ(a.k, b.k);
    ASSERT_EQ(a.N, b.N);
    ASSERT_EQ(a.b.coeffs, b.b.coeffs);
    ASSERT_EQ(a.b.N, b.b.N);
    for (auto i = 0; i < a.k; i++) {
        ASSERT_EQ(a.a[i].N, b.a[i].N);
        ASSERT_EQ(a.a[i].coeffs, b.a[i].coeffs);
    }
}

void verifyTrlweDft(const TrlweDft& a, const TrlweDft& b) {
    ASSERT_EQ(a.k, b.k);
    ASSERT_EQ(a.a_initialized, b.a_initialized);
    ASSERT_EQ(a.b.coeffs, b.b.coeffs);
    if (a.a_initialized) {
        for (auto i = 0; i < a.k; i++) {
            ASSERT_EQ(a.a[i].coeffs, b.a[i].coeffs);
        }
    }
}

void verifyTrlevDft(const TrlevDft& a, const TrlevDft& b) {
    ASSERT_EQ(a.l, b.l);
    for (auto i = 0; i < a.l; i++) {
        verifyTrlweDft(a.trlweDfts[i], b.trlweDfts[i]);
    }
}

void verifyTrgswMP(const TrgswMP& a, const TrgswMP& b) {
    ASSERT_EQ(a.l, b.l);
    ASSERT_EQ(a.k, b.k);
    ASSERT_EQ(a.isHalf, b.isHalf);
    for (auto l = 0; l < a.l; l++) {
        verifyTrlwe(a.cPrime[l],  b.cPrime[l]);
        if (a.isHalf) {
            continue;
        }
        for (auto i = 0; i < b.k; i++) {
            verifyTrlwe(a.c[l][i],  b.c[l][i]);
        }
    }
}

void verifyTrgswMPDft(const TrgswMPDft& a, const TrgswMPDft& b) {
    ASSERT_EQ(a.l, b.l);
    ASSERT_EQ(a.k, b.k);
    ASSERT_EQ(a.isHalf, b.isHalf);
    for (auto l = 0; l < a.l; l++) {
        verifyTrlweDft(a.cPrime[l],  b.cPrime[l]);
        if (a.isHalf) {
            continue;
        }
        for (auto i = 0; i < b.k; i++) {
            verifyTrlweDft(a.c[l][i],  b.c[l][i]);
        }
    }
}

void verifyBskMP(const BootstrappingKeyMP& a, const BootstrappingKeyMP& b) {
    ASSERT_EQ(a.n, b.n);
    ASSERT_EQ(a.group, b.group);
    ASSERT_EQ(a.isHalf, b.isHalf);
    for (auto i = 0; i < a.n; i++) {
        for (auto d = 0; d < a.bskDft[i].size(); d++) {
            verifyTrgswMPDft(a.bskDft[i][d], b.bskDft[i][d]);
        }
    }
}

void verifyBskMPOpt(const BootstrappingKeyMPOpt& a, const BootstrappingKeyMPOpt& b) {
    ASSERT_EQ(a.n, b.n);
    ASSERT_EQ(a.group, b.group);
    for (auto d = 0; d < a.bskFirst.size(); d++) {
        verifyTrlwe(a.bskFirst[d], b.bskFirst[d]);
    }
    for (auto i = 0; i < a.n - 1; i++) {
        for (auto d = 0; d < a.bskDft[i].size(); d++) {
            verifyTrgswMPDft(a.bskDft[i][d], b.bskDft[i][d]);
        }
    }
}

//void verifyBskMPLazy(const BootstrappingKeyMPLazy& a, const BootstrappingKeyMPLazy& b) {
//    ASSERT_EQ(a.n, b.n);
//    ASSERT_EQ(a.level, b.level);
//    ASSERT_EQ(a.group, b.group);
//    ASSERT_EQ(a.initialized, b.initialized);
//    for (auto d = 0; d < a.bskFirst.size(); d++) {
//        verifyTrlwe(a.bskFirst[d], b.bskFirst[d]);
//    }
//    for (auto i = 0; i < a.n - 1; i++) {
//        for (auto d = 0; d < a.bskTrim[i].size(); d++) {
//            verifyTrgswMPDft(a.bskTrim[i][d], b.bskTrim[i][d]);
//        }
//    }
//    for (auto i = 0; i < a.n - 1; i++) {
//        for (auto d = 0; d < a.bskDecompA[i].size(); d++) {
//            for (auto l1 = 0; l1 < a.bskDecompA[i][d].size(); l1++) {
//                for (auto l2 = 0; l2 < a.bskDecompA[i][d][l1].size(); l2++) {
//                    for (auto k = 0; k < a.bskDecompA[i][d][l1][l2].size(); k++) {
//                        ASSERT_EQ(a.bskDecompA[i][d][l1][l2][k].N, b.bskDecompA[i][d][l1][l2][k].N);
//                        ASSERT_EQ(a.bskDecompA[i][d][l1][l2][k].coeffs, b.bskDecompA[i][d][l1][l2][k].coeffs);
//                    }
//                }
//            }
//        }
//    }
//}

void verifyBskMPLazyAlt(const BootstrappingKeyMPLazyPipeAlt& a, const BootstrappingKeyMPLazyPipeAlt& b) {
    ASSERT_EQ(a.n, b.n);
    ASSERT_EQ(a.level, b.level);
    ASSERT_EQ(a.group, b.group);
    verifyTrlevDft(a.bskFirstLev, b.bskFirstLev);
    verifyTrlevDft(a.s2Dft, b.s2Dft);
    for (auto i = 0; i < a.n - 1; i++) {
        verifyTrgswMP(a.bskPrime[i], b.bskPrime[i]);
    }
}

void verifyBskWWL24(const BootstrappingKeyWWL24& a, const BootstrappingKeyWWL24& b) {
    ASSERT_EQ(a.n, b.n);
    ASSERT_EQ(a.group, b.group);
    verifyTrlevDft(a.s2Dft, b.s2Dft);
    for (auto i = 0; i < a.n; i++) {
        verifyTrgswMPDft(a.bskDft[i], b.bskDft[i]);
    }
}

TEST(SERIALIZATION, POLYNOMIAL) {
    YatfheParameters param {};
    param.N = 128;
    initYatfhe(param);

    TorusPolynomial polyT32A{param.N};
    NttPolynomial ntt{param.N};
    DecompPolynomial dp{param.N};

    for (auto i = 0; i < param.N; i++) {
        polyT32A.coeffs[i] = modSwitchToTorus32(intModP(i, param.torusBase), param.torusBase);
        dp.coeffs[i] = longModP(i*(1<<12), 1<<16);
    }
    NttHexl::applyNtt(ntt, polyT32A);

    printArray(polyT32A.coeffs, "polyT32A");
    printArray(ntt.coeffs, "ntt");
    printArray(dp.coeffs, "dp");


    std::ostringstream oss(std::ios::binary);
    vector<char> buffer;
    serialize(polyT32A, oss);
    const std::string &data = oss.str();
    buffer.insert(buffer.end(), data.begin(), data.end());
    oss.str("");
    oss.clear();
    serialize(ntt, oss);
    const std::string &data1 = oss.str();
    buffer.insert(buffer.end(), data1.begin(), data1.end());
    oss.str("");
    oss.clear();
    serialize(dp, oss);
    const std::string &data2 = oss.str();
    buffer.insert(buffer.end(), data2.begin(), data2.end());

    std::ofstream outFile("SERIALIZATION_TEST_POLY.bin", std::ios::binary | std::ios::trunc);
    outFile.write(buffer.data(), buffer.size());
    buffer.clear();
    outFile.close();

    TorusPolynomial read;
    NttPolynomial readNtt;
    DecompPolynomial readDp;
    std::ifstream inFile("SERIALIZATION_TEST_POLY.bin", std::ios::binary);
    deserialize(read, inFile);
    deserialize(readNtt, inFile);
    deserialize(readDp, inFile);
    inFile.close();
    printArray(read.coeffs, "read");
    printArray(readNtt.coeffs, "readNtt");
    printArray(readDp.coeffs, "readDp");

    ASSERT_EQ(read.coeffs, polyT32A.coeffs);
    ASSERT_EQ(readNtt.coeffs, ntt.coeffs);
    ASSERT_EQ(readDp.coeffs, dp.coeffs);

    printBanner("SERIALIZATION.POLYNOMIAL");
}

TEST(SERIALIZATION, TRLWE) {
    YatfheParameters param{};
    param.N = 64;
    initYatfhe(param);

    TrlweKey trlweKey{param};
    Trlwe trlwe{param};
    TrlweDft trlweDft{param};
    genTrlweKey(trlweKey);

    TorusPolynomial torus{param.N};
    for (auto i = 0; i < param.N; i++) {
        torus.coeffs[i] =  modSwitchToTorus32(intModP(i, param.torusBase), param.torusBase);;
    }

    symEncTrlweMultiSampleNtt(trlwe, trlweDft, trlweKey, torus.coeffs);
    printTrlweAB(trlwe, "trlwe");
    printTrlweDftAB(trlweDft, "trlweDft");

    std::ostringstream oss(std::ios::binary);
    vector<char> buffer;
    serialize(trlwe, oss);
    const std::string &data = oss.str();
    buffer.insert(buffer.end(), data.begin(), data.end());
    oss.str("");
    oss.clear();
    serialize(trlweDft, oss);
    const std::string &data1 = oss.str();
    buffer.insert(buffer.end(), data1.begin(), data1.end());

    std::ofstream outFile("SERIALIZATION_TEST_TRLWE.bin", std::ios::binary | std::ios::trunc);
    outFile.write(buffer.data(), buffer.size());
    buffer.clear();
    outFile.close();

    Trlwe readTrlwe;
    TrlweDft readTrlweDft;
    std::ifstream inFile("SERIALIZATION_TEST_TRLWE.bin", std::ios::binary);
    deserialize(readTrlwe, inFile);
    deserialize(readTrlweDft, inFile);
    inFile.close();

    printTrlweAB(readTrlwe, "readTrlwe");
    printTrlweDftAB(readTrlweDft, "readTrlweDft");

    verifyTrlwe(trlwe, readTrlwe);
    verifyTrlweDft(trlweDft, readTrlweDft);

    printBanner("SERIALIZATION.TRLWE");

}

TEST(SERIALIZATION, TRGSW) {
    YatfheParameters param{};
    param.N = 64;
    initYatfhe(param);

    // key gen
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

    // trgsw enc
    TrgswMPDft trgswDft{param, param.lApprox};
    Integer mu1 = genIntUniformDist(0, 3);
    encryptTrgswMPNtt(trgswDft, mu1, trgswKey, 0, param);

    std::ostringstream oss(std::ios::binary);
    vector<char> buffer;
    serialize(trgswDft, oss);
    const std::string &data = oss.str();
    buffer.insert(buffer.end(), data.begin(), data.end());
    std::ofstream outFile("SERIALIZATION_TEST_TRGSW.bin", std::ios::binary | std::ios::trunc);
    outFile.write(buffer.data(), buffer.size());
    buffer.clear();
    outFile.close();

    TrgswMPDft readTrgswDft;
    COUNT_TIME("inFile", std::ifstream inFile("SERIALIZATION_TEST_TRGSW.bin", std::ios::binary);)
    COUNT_TIME("deserialize", deserialize(readTrgswDft, inFile);)
    COUNT_TIME("close", inFile.close();)
    verifyTrgswMPDft(readTrgswDft, trgswDft);

    printBanner("SERIALIZATION.TRGSW");

}

TEST(SERIALIZATION, BSKMP) {
    YatfheParameters param{};
    param.N = 1024;
    initYatfhe(param);

    // client side
    // key gen
    TlweKey tlweKey{param.n};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TrlevDft s2Dft{param, param.l};
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
    serializeBskMP(bskMP, "SERIALIZATION_TEST_BSKMP.bin");

    BootstrappingKeyMP readBskMP;
    deserializeBskMP(readBskMP, "SERIALIZATION_TEST_BSKMP.bin", param.n);

    verifyBskMP(bskMP, readBskMP);

    printBanner("SERIALIZATION.BSKMP");
}

TEST(SERIALIZATION, BSKMPOPT) {
    YatfheParameters param{};
    param.N = 1024;
    initYatfhe(param);

    // client side
    // key gen
    TlweKey tlweKey{param.n};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TrlevDft s2Dft{param, param.l};
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    serializeBskMPOpt(bskMPOpt, "SERIALIZATION_TEST_BSKMPOPT.bin");

    BootstrappingKeyMPOpt readBskMPOpt;
    deserializeBskMPOpt(readBskMPOpt, "SERIALIZATION_TEST_BSKMPOPT.bin", param.n);

    verifyBskMPOpt(bskMPOpt, readBskMPOpt);

    printBanner("SERIALIZATION.BSKMPOPT");
}

//TEST(SERIALIZATION, BSKMPLAZY) {
//    YatfheParameters param{};
//    initYatfhe(param);
//
//    // client side
//    // key gen
//    TlweKey tlweKey{param};
//    TrgswKey trgswKey{param};
//    TrlweKey& trlweKey = trgswKey.trlweKey;
//    TlweKeySwitchingKey ksKey{param};
//    genTlweKey(tlweKey);
//    genTrlweKey(trlweKey);
//
//    TrlevDft s2Dft{param, param.l};
//    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
//    TorusPolynomial v {param.N};
//    generateTestPolynomial(v, param.torusBase, 2 * param.N);
//
//    BootstrappingKeyMPLazy bskMPLazy{param, param.lApprox, true, true};
//    genBootstrappingKeyMPLazy(bskMPLazy, trgswKey, tlweKey, v, param);
//
//    COUNT_TIME("serializeBskMPLazy", serializeBskMPLazy(bskMPLazy, "SERIALIZATION_TEST_BSKMPLAZY.bin");)
//
//    BootstrappingKeyMPLazy readBskMPLazy;
//    COUNT_TIME("deserializeBskMPLazy", deserializeBskMPLazy(readBskMPLazy, "SERIALIZATION_TEST_BSKMPLAZY.bin", param.n);)
//
//    verifyBskMPLazy(bskMPLazy, readBskMPLazy);
//
//    printBanner("SERIALIZATION.BSKMPLAZY");
//}

TEST(SERIALIZATION, BSKMPLAZY_ALT) {
    YatfheParameters param{};
    initYatfhe(param);

    // client side
    // key gen
    TlweKey tlweKey{param};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
    symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trlweKey, trlweKey.s, true, param);
    genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, param);

    COUNT_TIME("serializeBskLazyPipeAlt", serializeBskLazyPipeAlt(bskMPLazyPipeAlt, "SERIALIZATION_TEST_BSKMPLAZY_ALT.bin");)

    BootstrappingKeyMPLazyPipeAlt readBskMPLazyAlt;
    COUNT_TIME("deserializeBskLazyPipeAlt", deserializeBskLazyPipeAlt(readBskMPLazyAlt, "SERIALIZATION_TEST_BSKMPLAZY_ALT.bin", param.n);)

    verifyBskMPLazyAlt(bskMPLazyPipeAlt, readBskMPLazyAlt);

    printBanner("SERIALIZATION.BSKMPLAZY_ALT");
}

TEST(SERIALIZATION, WWL24) {
    YatfheParameters param{};
    initYatfhe(param);

    // client side
    // key gen
    TlweKey tlweKey{param};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyWWL24 bskWWL24{param, param.lApprox};
    COUNT_TIME("genBootstrappingKeyWWL24", genBootstrappingKeyWWL24(bskWWL24, trgswKey, tlweKey, param);)

    COUNT_TIME("serializeBskWWL24", serializeBskWWL24(bskWWL24, "SERIALIZATION_TEST_BSKWWL24.bin");)

    BootstrappingKeyWWL24 readBskWWL24;
    COUNT_TIME("deserializeBskWWL24", deserializeBskWWL24(readBskWWL24, "SERIALIZATION_TEST_BSKWWL24.bin", param.n);)

    verifyBskWWL24(bskWWL24, readBskWWL24);

    printBanner("SERIALIZATION.WWL24");
}