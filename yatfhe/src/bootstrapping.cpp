//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/ntt24.h"

void trgswFunctionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, 0); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    tlweKeySwitch(out, ksk, tmp, param);
}

void trgswFunctionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotateNtt(accum, bsk, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, 0); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    tlweKeySwitch(out, ksk, tmp, param);
}

//todo
void trgswFunctionalBootstrappingCRT(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(tv, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    trlweMCRTDecomp(accCRT, tv, param);
    blindRotateMCRT(accCRT, bskCRT, inputModN2, param);
    trlweCRTRecomp(acc, accCRT, param);
    extractTlweFromTrlwe(tmp, acc, 0); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    tlweKeySwitch(out, ksk, tmp, param);
}

/**
 * Multiply the accumulator by X^sum(bara_i * s_i)
 * */
void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = Trlwe{param.k, param.N};
        controlMux(temp, accum, input.a[i], bsk.bsk[i], param);
        swap(accum, temp); // assign the previous result to accumulator
    }
}

void blindRotateNtt(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = Trlwe{param.k, param.N};
        controlMuxNtt(temp, accum, input.a[i], bsk.bskDft[i], param);
        swap(accum, temp); // assign the previous result to accumulator
    }
}

void blindRotateMCRT(std::vector<Trlwe8>& accum, const BootstrappingKeyCRT& bskCRT, const ScaledTlwe& input, const YatfheParameters& param) {
    std::vector<Trlwe8> temp(param.d, Trlwe8{param.k, param.N});
    for (size_t i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        for (size_t d = 0; d < param.d; d++) {
            temp[d] = Trlwe8{param.k, param.N};
        }
        controlMuxCRT(temp, accum, input.a[i], bskCRT.bskCRT[i], param);
        swap(accum, temp);
    }
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMux(Trlwe& res, const Trlwe& input, const int aBarI, const Trgsw& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    trlweRotateMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    trgswExternalProduct(res, bskI, tmp, param); // res *= bskI
    trlweAccumulateT32(res, input); // res += input
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMuxNtt(Trlwe& res, const Trlwe& input, const int aBarI, const TrgswDft& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    trlweRotateMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    trgswExternalProductNtt(res, bskI, tmp, param); // res *= bskI
    trlweAccumulateT32(res, input); // res += input
}

void controlMuxCRT(std::vector<Trlwe8>& res, const std::vector<Trlwe8>& inputs, const int aBarI, const std::vector<TrgswDft24>& bskCRT, const YatfheParameters& param) {
    std::vector<Trlwe8> tmp(param.d, Trlwe8{param.k, param.N});

    for (size_t i = 0; i < param.d; i++) {
        trlweRotateMinusOne8(tmp[i], inputs[i], aBarI, param.qd[i]); // res = c1 - c0 = X^aBarI * input - input
    }

    trgswExternalProductCRTNTT(res, bskCRT, tmp, param); // res *= bskI

    for (size_t i = 0; i < param.d; i++) {
        trlweAccumulateModP(res[i], inputs[i], param.qd[i]); // res += input
    }
}

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey) {
    if (bsk.unfolding == 1) {
        bootstrappingKeyGenWoUnfolding(bsk, param, trgswKey, tlweKey);
    }
}

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey) {
    for (auto i = 0; i < bsk.n; i++) {
        trgswEncryptNtt(bsk.bsk[i], bsk.bskDft[i], param, trgswKey, tlweKey.s[i]);
    }
}

void bootstrappingKeyMCRTDecomp(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param) {
    for (size_t i = 0; i < param.n; i++) {
        trgswMCRTDecomp(bskCRT.bsk8[i], bsk.bsk[i], param);
        auto& bsk8D = bskCRT.bsk8[i];
        auto& bskNttD = bskCRT.bskCRT[i];
        for (size_t d = 0; d < param.d; d++) {
            auto& bsk8 = bsk8D[d];
            auto& bskNtt = bskNttD[d];
            for (size_t l = 0; l < param.dh; l++) {
                for (size_t k1 = 0; k1 < param.k + 1; k1++) {
                    auto& nttOut = bskNtt.trlweDftSamples[l][k1];
                    auto& nttIn = bsk8.trlweSamples[l][k1];
                    applyNttForAB24(nttOut, nttIn);
                }
            }
        }
    }
}














