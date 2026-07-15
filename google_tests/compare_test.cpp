//
// Created by xsong93 on 3/26/26.
//

#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/ntt_hexl.h"

TEST(COMPARE, MP21_NTT) {
    YatfheParameters param {};
    constexpr int kPlainMod = 1 << 4;
    param.torusBase = kPlainMod;
    param.driftPhase = param.N / param.torusBase / 2;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey {param};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyMP bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    genBootstrappingKeyMP(bsKey, trgswKey, tlweKey, param);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    int threshold = kPlainMod / 4;
    cout << "threshold: " << threshold << endl;
    TorusPolynomial v {param.N};
    generateTestPolynomialLt1(v, threshold);

    for (int plain = 0; plain < kPlainMod / 2; plain += 1) {
        Torus mu = modSwitchToTorusGeneral(plain, kPlainMod, LWE_Q);

        Tlwe input {param.n};
        Tlwe output {param.n};
        symEncTlwe(input, mu, tlweKey);

        auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
        cout << "decPre: " << decPre << endl;

        ScaledTlwe inputModN2 {param.N * 2, param.n};
        Trlwe accum {param.k, param.N};
        Tlwe tmp {ksKey.nCurrKey};

        rescaleTlweToNewMod(inputModN2, input);
        genNoiselessTrlweSample(accum, v, inputModN2);
        blindRotateMP21Ntt(accum, bsKey.bskDft, inputModN2, param);
        extractTlweFromTrlwe(tmp, accum, param.driftPhase);
        switchKeyForTlwe(output, ksKey, tmp, param);

        auto decAft = symDecTlweToInt(output, tlweKey, kPlainMod);
        const int expected = plain < threshold ? 1 : 0;
        ASSERT_EQ(expected, decAft);
        cout << "plain: " << plain << ", decAft: " << decAft << endl;
    }

    printBanner("COMPARE.MP21_NTT");
}

TEST(COMPARE, COMP_ZERO) {
    YatfheParameters param {};
    constexpr int kPlainMod = 1 << 4;
    param.torusBase = kPlainMod;
    param.driftPhase = param.N / param.torusBase / 2;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey {param};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyMP bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    genBootstrappingKeyMP(bsKey, trgswKey, tlweKey, param);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    int threshold = 0;
    Integer value = 2;
    cout << "threshold: " << threshold << endl;
    TorusPolynomial v {param.N};
    generateTestPolynomialCompWithValue(v, threshold, value);

    for (int plain = -kPlainMod / 2 + 1; plain < kPlainMod / 2; plain += 1) {
        Torus mu = modSwitchToTorusGeneral(plain, kPlainMod, LWE_Q);

        Tlwe input {param.n};
        Tlwe output {param.n};
        symEncTlwe(input, mu, tlweKey);

        auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
        cout << "decPre: " << decPre << endl;

        ScaledTlwe inputModN2 {param.N * 2, param.n};
        Trlwe accum {param.k, param.N};
        Tlwe tmp {ksKey.nCurrKey};

        rescaleTlweToNewMod(inputModN2, input);
        genNoiselessTrlweSample(accum, v, inputModN2);
        blindRotateMP21Ntt(accum, bsKey.bskDft, inputModN2, param);
        extractTlweFromTrlwe(tmp, accum, param.driftPhase);
        switchKeyForTlwe(output, ksKey, tmp, param);

        auto decAft = symDecTlweToInt(output, tlweKey, kPlainMod);
        const int expected = plain == threshold ? value : 0;
        cout << "plain: " << plain << ", decAft: " << decAft << endl;
        ASSERT_EQ(expected, decAft);
    }

    printBanner("COMPARE.COMP_ZERO");
}