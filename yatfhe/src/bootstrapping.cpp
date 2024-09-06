//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/ntt.h"
#include "yatfhe/crt.h"
#include "yatfhe/ntt14.h"

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

/**
 * Multiply the accumulator by X^sum(bara_i * s_i)
 * */
void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe temp {param.k, param.N};
        controlMux(temp, accum, input.a[i], bsk.bsk[i], param);
        swap(accum, temp); // assign the previous result to accumulator
    }
}

void blindRotateNtt(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe temp {param.k, param.N};
        controlMuxNtt(temp, accum, input.a[i], bsk.bskDft[i], param);
        swap(accum, temp); // assign the previous result to accumulator
    }
}

//todo
void blindRotateCRT(Trlwe8& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe8 temp {param.k, param.N};
//        controlMuxCRT(temp, accum, input.a[i], bsk.bskDft[i], param);
        swap(accum, temp); // assign the previous result to accumulator
    }
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMux(Trlwe& res, const Trlwe& input, const int aBarI, const Trgsw& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    trlweRotateMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    trgswExternalProduct(res, bskI, tmp, param); // res *= bskI
    trlweAccumulate(res, input); // res += input
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMuxNtt(Trlwe& res, const Trlwe& input, const int aBarI, const TrgswDft& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    trlweRotateMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    trgswExternalProductNtt(res, bskI, tmp, param); // res *= bskI
    trlweAccumulate(res, input); // res += input
}

// todo
void controlMuxCRT(Trlwe8& res, const std::vector<Trlwe8>& inputs, const int aBarI, const TrgswDft& bskI, const YatfheParameters& param) {
    std::vector<Trlwe8> tmp (param.d, Trlwe8(param.k, param.N));
    std::vector<Trlwe8> tmpD (param.dh, Trlwe8(param.k, param.N));
    std::vector<std::vector<Trlwe8>> tmpDB (param.d, std::vector<Trlwe8>(param.dh, Trlwe8(param.k, param.N)));
    std::vector<std::vector<TrlweDft14>> tmpDBNtt (param.d, std::vector<TrlweDft14>(param.dh, TrlweDft14(param.k, param.N)));

    for (size_t i = 0; i < param.d; i++) {
        trlweRotateMinusOne8(tmp[i], inputs[i], aBarI, param.qd[i]); // res = c1 - c0 = X^aBarI * input - input
    }
    syncGadgetDecomp(tmpD, tmp, param);
    broadcastCRT(tmpDB, tmpD, param);

    // ntt
    for (size_t i1 = 0; i1 < param.d; i1++) {
        for (size_t i2 = 0; i2 < param.dh; i2++) {
            for (size_t k = 0; k < param.k; k++) {
                applyNtt14Poly8(tmpDBNtt[i1][i2].a[k], tmpDB[i1][i2].a[k]);
            }
            applyNtt14Poly8(tmpDBNtt[i1][i2].b, tmpDB[i1][i2].b);
        }
    }
//    trgswExternalProductCRT(res, bskI, tmp, param); // res *= bskI
//    trlweAccumulate(res, input); // res += input
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














