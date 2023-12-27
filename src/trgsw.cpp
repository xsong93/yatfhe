//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "trgsw.h"
#include "trlwe.h"

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const int l, const int bgBit) {
    trgswKey.trlweKey = &trlweKey;
    trgswKey.l = l;
    trgswKey.bgBit = bgBit;
}

void initTrgswSample(Trgsw& trgsw, const int l, const int bgBit, const int k, const int N) {
    for (int i = 0; i < l * (k + 1); i++) {
        initTrlweSample(trgsw.trlweSamples[i], k, N);
    }
}

void initTrgswDftSample(TrgswDft& trgswDft, const int l, const int bgBit, const int k, const int N) {
    for (size_t i = 0; i < l * (k + 1); i++) {
        initTrlweDftSample(trgswDft.trlweDftSamples[i], k, N);
    }
}

void deleteTrgswKey(TrgswKey& trgswKey) {
    trgswKey.trlweKey = nullptr;
}






