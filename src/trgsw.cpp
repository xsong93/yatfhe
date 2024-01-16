//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "trgsw.h"
#include "trlwe.h"

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const int l, const int bgBit) {
    trgswKey.trlweKey = trlweKey;
    trgswKey.l = l;
    trgswKey.bgBit = bgBit;
}

void initTrgswSample(Trgsw& trgsw, const TrgswKey& trgswKey) {
    const int l = trgswKey.l;
    const int k = trgswKey.trlweKey.k;
    const int N = trgswKey.trlweKey.s[0].N;
    const int kpl = l * (k + 1);
    trgsw.trlweSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweSample(trgsw.trlweSamples[i], k, N);
    }
}

void initTrgswDftSample(TrgswDft& trgswDft, const TrgswKey& trgswKey) {
    const int l = trgswKey.l;
    const int k = trgswKey.trlweKey.k;
    const int N = trgswKey.trlweKey.s[0].N;
    const int kpl = l * (k + 1);
    trgswDft.trlweDftSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        initTrlweDftSample(trgswDft.trlweDftSamples[i], k, N);
    }
}

void allocNewTrgswDftSample(TrgswDft& trgswDftSample, const int l, const int Bg_bit, const int k, const int N) {
    const int kpl = l *  (k + 1);
    trgswDftSample.trlweDftSamples.resize(kpl);
    for (int i = 0; i < kpl; i++) {
        trgswDftSample.trlweDftSamples[i] = trlwe_alloc_new_DFT_sample(k, N);
    }
}

void deleteTrgswKey(TrgswKey& trgswKey) {
//    trgswKey.trlweKey = nullptr;
}






