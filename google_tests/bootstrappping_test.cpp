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
#include "yautil/time_counter.h"
#include "yatfhe/blind_rotate.h"

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
    TlweKey tlweKey {param};
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
    printBanner("BOOTSTRAPPING.GROUP2_KEYGEN");
}

TEST(BOOTSTRAPPING, MP21) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey {param.n};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyMP bsKey {param, param.lApprox};
    TlweKeySwitchingKey ksKey {param};
    COUNT_TIME("genTlweKey", genTlweKey(tlweKey);)
    COUNT_TIME("genTrlweKey", genTrlweKey(trlweKey);)
    COUNT_TIME("genBootstrappingKey", genBootstrappingKeyMP(bsKey, trgswKey, tlweKey, param);)
    TlweKey tlweKsKey = tlweKey;
    COUNT_TIME("genTlweKeySwitchingKey", genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);)

    Integer plain = 3;
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorusGeneral(mu, param.torusBase, LWE_Q) << endl;
    auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Trlwe accumScaled {param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    COUNT_TIME("rescaleTlweToNewMod", rescaleTlweToNewMod(inputModN2, input);) // rescale to mod 2N
    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(accum, v, inputModN2);) // accum = (X^-b) * (0,...,0,v)
    COUNT_TIME("blindRotateMP21Ntt", blindRotateMP21Ntt(accum, bsKey.bskDft, inputModN2, param);)
    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, accum, param.driftPhase);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    COUNT_TIME("switchKeyForTlwe", switchKeyForTlwe(output, ksKey, tmp, param);)

    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    ASSERT_EQ(decPre, decAft);
    cout << "err:" << calTlweError(output, tlweKey, plain) << endl;

    printBanner("BOOTSTRAPPING.MP21");
}

TEST(BOOTSTRAPPING, MP21_FR) {
    YatfheParameters param {};
    param.torusBase = 8;
    // param.setRadixBits(8);
    // param.l = 4;
    // param.lApprox = 3;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey {param};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyMP bskMP{param, param.lApprox};
    TlweKeySwitchingKey ksKey {param};
    COUNT_TIME("genTlweKey", genTlweKey(tlweKey);)
    COUNT_TIME("genTrlweKey", genTrlweKey(trlweKey);)
    COUNT_TIME("genBootstrappingKey", genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);)
    TlweKey tlweKsKey = tlweKey;
    COUNT_TIME("genTlweKeySwitchingKey", genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);)

    int plain = -3;
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    TorusPolynomial v {param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorusGeneral(mu, param.torusBase, LWE_Q) << endl;
    auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

    ScaledTlwe inputModN2 {param.N, param.n};
    Trlwe accum {param};
    Trlwe accumScaled {param};
    Tlwe tmp {ksKey.nCurrKey};
    COUNT_TIME("rescaleTlweToNewMod", rescaleTlweToNewMod(inputModN2, input);) // rescale to mod 2N
    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(accum, v, inputModN2);) // accum = (X^-b) * (0,...,0,v)
    COUNT_TIME("blindRotateJP22Ntt", blindRotateJP22Ntt(accum, bskMP, inputModN2, param);)
    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, accum, 0);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    COUNT_TIME("switchKeyForTlwe", switchKeyForTlwe(output, ksKey, tmp, param);)

    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    ASSERT_EQ(decPre, decAft);
    cout << "err:" << calTlweError(output, tlweKey, plain) << endl;

    printBanner("BOOTSTRAPPING.MP21_FR");
}

TEST(BOOTSTRAPPING, MCRT) {
    YatfheParameters param{};
    param.q = Q_CRT;
    param.torusBits = 32;
    param.torusBase = 8;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyCRT bsKeyCRT{param};
    TlweKeySwitchingKey ksKey{param};
    COUNT_TIME("genTlweKey", genTlweKey(tlweKey);)
    COUNT_TIME("genTrlweKey", genTrlweKey(trlweKey);)
    COUNT_TIME("genBootstrappingKeyApproxCrt", genBootstrappingKeyApproxCrt(bsKeyCRT, trgswKey, tlweKey, param);)
    TlweKey tlweKsKey = tlweKey;
    COUNT_TIME("genTlweKeySwitchingKey", genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);)

    Integer plain = 3;
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    TorusPolynomial v{param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    Tlwe input{param.n};
    Tlwe output{param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorusGeneral(mu, param.torusBase, LWE_Q) << endl;
    auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksKey.nCurrKey};
//    COUNT_TIME("rescaleTlweToNewMod", rescaleTlweToNewMod(inputModN2, input);) // rescale to mod 2N
//    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(tv, v, inputModN2);) // tv = (X^-b) * (0,...,0,v)
//    COUNT_TIME("decompTrlweMcrt", decompTrlweMcrt(accCRT, tv, param);)
//    COUNT_TIME("blindRotateApproxCRTNtt", blindRotateApproxCRTNtt(accCRT, bsKeyCRT.bskCRT, inputModN2, param);)
//    COUNT_TIME("trlweMcrtToCrt", trlweMcrtToCrt(accCRT, param);)
//    COUNT_TIME("recompTrlweCrt", recompTrlweCrt(acc, accCRT, param);)
//    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, acc, param.driftPhase);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
//    COUNT_TIME("switchKeyForTlwe", switchKeyForTlwe(output, ksKey, tmp, param);)

    COUNT_TIME("functionalBootstrappingCrt", functionalBootstrappingCrt(output, input, bsKeyCRT, ksKey, v, param);)

    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    ASSERT_EQ(decPre, decAft);
    cout << "err:" << calTlweError(output, tlweKey, plain) << endl;

    printBanner("BOOTSTRAPPING.MCRT");
}