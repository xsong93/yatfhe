//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "trgsw.h"
#include "trlwe.h"

TrgswKey* trgswInitKey(TrlweKey* trlweKey, const int l, const int bgBit) {
    TrgswKey* res{new TrgswKey};
    res->trlweKey = trlweKey;
    res->l = l;
    res->bgBit = bgBit;
    return res;
}

void initTrgswSample(Trgsw* trgsw, const int l, const int bgBit, const int k, const int N) {
    trgsw->samples = new Trlwe;
    for (int i = 0; i < l * (k + 1); i++) {
//        res->samples[i] = *initTrlweSample(k, N);
        initTrlweSample(&trgsw->samples[i], k, N);
    }
//    trgsw->bgBit = bgBit;
//    trgsw->l = l;
}
void initTrgswDftSample(TrgswDft* trgswDft, const int n, const int l, const int bgBit, const int k, const int N) {
    for (int i = 0; i < n; i++) {
//        res[i] = *initTrgswDftSample(l, bgBit, k, N);
        initTrgswDftSample(&trgswDft[i], l, bgBit, k, N);
    }
}

void initTrgswDftSample(TrgswDft* trgsw, const int l, const int bgBit, const int k, const int N) {
    trgsw->samples = new TrlweDft;
    for (size_t i = 0; i < l * (k + 1); i++) {
//        trgsw->samples[i] = *initTrlweDftSample(k, N);
        initTrlweDftSample(&trgsw->samples[i], k, N);
    }
//    trgsw->bgBit = bgBit;
//    trgsw->l = l;
}

//TrgswDft* initTrgswDftSample(int l, int bgBit, int k, int N) {
//    TrgswDft* res{new TrgswDft};
//    res->samples = new TrlweDft;
//    for (size_t i = 0; i < l * (k + 1); i++) {
//        res->samples[i] = *initTrlweDftSample(k, N);
//    }
//    res->bgBit = bgBit;
//    res->l = l;
//    return res;
//}