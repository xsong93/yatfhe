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
    trgsw.trlweSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweSample(trgsw.trlweSamples[i], k, N);
    }
}

void initTrgswDftSample(TrgswDft& trgswDftSample, const YatfheParameters& yatfheParameters) {
    const int l = yatfheParameters.l;
    const int k = yatfheParameters.k;
    const int N = yatfheParameters.N;
    const int kpl = l *  (k + 1);
    trgswDftSample.trlweDftSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweDftSample(trgswDftSample.trlweDftSamples[i], k, N);
    }
}

void deleteTrgswKey(TrgswKey& trgswKey) {
//    trgswKey.trlweKey = nullptr;
}






