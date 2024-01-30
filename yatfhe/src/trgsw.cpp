//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe_parameters.h"
#include "trgsw.h"
#include "trlwe.h"

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const YatfheParameters& param) {
    trgswKey.trlweKey = trlweKey;
    trgswKey.l = param.l;
    trgswKey.bgBit = param.bgBit;
}

void initTrgswSample(Trgsw& trgsw, const YatfheParameters& param) {
    const int l = param.l;
    const int k = param.k;
    const int N = param.N;
    const int kpl = l * (k + 1);
//    trgsw.trlweSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweSample(trgsw.trlweSamples[i], k, N);
    }
}

void initTrgswDftSample(TrgswDft& trgswDftSample, const YatfheParameters& param) {
    const int l = param.l;
    const int k = param.k;
    const int N = param.N;
    const int kpl = l *  (k + 1);
//    trgswDftSample.trlweDftSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweDftSample(trgswDftSample.trlweDftSamples[i], k, N);
    }
}

void genNoiselessTrgswSample(Trgsw &trgswSample, Torus msg, const YatfheParameters& param) {
    const int l = param.l;
    const int bgBit = param.bgBit;
    const int k = param.k;
    const int N = param.N;
    const int kpl = (k + 1) * l;
    for (int i = 0; i < l; i++) {
        const uint64_t h = 1UL << (sizeof(Torus) * 8 - (i + 1) * bgBit);
        for (int j = 0; j < k; j++) {
            trgswSample.trlweSamples[j * l + i].a[j].coeffs[0] += msg * h;
        }
        trgswSample.trlweSamples[k * l + i].b.coeffs[0] += msg * h;
    }
}

void deleteTrgswKey(TrgswKey& trgswKey) {
//    trgswKey.trlweKey = nullptr;
}






