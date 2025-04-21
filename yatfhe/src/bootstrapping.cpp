//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/key_patterns.h"

void functionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk.bsk, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotateNtt(accum, bsk.bskDft, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingCrt(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input);// rescale to mod 2N
    genNoiselessTrlweSample(tv, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    decompTrlweMcrt(accCRT, tv, param);
    blindRotateApproxCRTNtt(accCRT, bskCRT.bskCRT, inputModN2, param);
    trlweMcrtToCrt(accCRT, param);
    recompTrlweCrt(acc, accCRT, param);
    extractTlweFromTrlwe(tmp, acc, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void genBootstrappingKeyGroup2(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    int j = 0;
    auto batchSize = 1 << param.group;
    for (auto i = 0; i < tlweKey.n; i = i + 2) {
        auto s1 = tlweKey.s[i];
        auto s2 = tlweKey.s[i + 1];
        int combined = (s1 << 1) | s2;

#pragma unroll(4)
        for (int k = 0; k < batchSize; ++k) {
            encryptTrgswNtt(bsk.bsk[j + k], bsk.bskDft[j + k], KEY_PATTERNS2[combined][k], trgswKey, 0, param);
        }
        j += batchSize;
    }
}

void genBootstrappingKeyNormal(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    for (auto i = 0; i < bsk.n; i++) {
        encryptTrgswNtt(bsk.bsk[i], bsk.bskDft[i], tlweKey.s[i], trgswKey, 0, param);
    }
}

void genBootstrappingKey(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    switch(bsk.group) {
        case 2:
            genBootstrappingKeyGroup2(bsk, trgswKey, tlweKey, param);
            return;
        default:
            genBootstrappingKeyNormal(bsk, trgswKey, tlweKey, param);
    }
}

void genBootstrappingKeyApproxCrt(BootstrappingKeyCRT& bskCRT, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    BootstrappingKey bsk{param};
    for (auto i = 0; i < bsk.n; i++) {
        encryptTrgswApproxCRT(bsk.bsk[i], param, trgswKey, tlweKey.s[i]);
    }
    decompBootstrappingKeyMcrt(bskCRT, bsk, param);
}

void decompBootstrappingKeyMcrt(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param) {
    for (size_t i = 0; i < param.n; i++) {
        decompTrgswMcrt(bskCRT.bsk8[i], bsk.bsk[i], param);
        auto& bsk8D = bskCRT.bsk8[i];
        auto& bskNttD = bskCRT.bskCRT[i];
        for (size_t d = 0; d < param.d; d++) {
            auto& bsk8 = bsk8D[d];
            auto& bskNtt = bskNttD[d];
            NttNative24::applyNttForRgsw(bskNtt, bsk8);
        }
    }
}