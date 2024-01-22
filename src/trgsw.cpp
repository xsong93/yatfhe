//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe_parameters.h"
#include "trgsw.h"
#include "trlwe.h"

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const YatfheParameters& yatfheParameters) {
    trgswKey.trlweKey = trlweKey;
    trgswKey.l = yatfheParameters.l;
    trgswKey.bgBit = yatfheParameters.bgBit;
}

void initTrgswSample(Trgsw& trgsw, const YatfheParameters& yatfheParameters) {
    const int l = yatfheParameters.l;
    const int k = yatfheParameters.k;
    const int N = yatfheParameters.N;
    const int kpl = l * (k + 1);
//    trgsw.trlweSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweSample(trgsw.trlweSamples[i], k, N);
    }
}

void initTrgswDftSample(TrgswDft& trgswDftSample, const YatfheParameters& yatfheParameters) {
    const int l = yatfheParameters.l;
    const int k = yatfheParameters.k;
    const int N = yatfheParameters.N;
    const int kpl = l *  (k + 1);
//    trgswDftSample.trlweDftSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweDftSample(trgswDftSample.trlweDftSamples[i], k, N);
    }
}

void genNoiselessTrgswSample(Trgsw &trgswSample, Torus msg, const YatfheParameters& parameters) {
    const int l = parameters.l;
    const int bgBit = parameters.bgBit;
    const int k = parameters.k;
    const int N = parameters.N;
    const int kpl = (k + 1) * l;
    for (int i = 0; i < l; i++) {
        const uint64_t h = 1UL << (sizeof(Torus)*8 - (i + 1) * bgBit);
        for (int j = 0; j < k; j++) {
            trgswSample.trlweSamples[j * l + i].a[j].coeffs[0] += msg * h;
        }
        trgswSample.trlweSamples[k * l + i].b.coeffs[0] += msg * h;
    }
}

void deleteTrgswKey(TrgswKey& trgswKey) {
//    trgswKey.trlweKey = nullptr;
}






