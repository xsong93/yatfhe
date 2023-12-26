//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include "trlwe.h"

struct Trgsw {
    Trlwe* samples; // l * (k + 1)
//    int l;
//    int bgBit;
};

struct TrgswKey {
    TrlweKey* trlweKey;
    int l;
    int bgBit;
};

struct TrgswDft {
    TrlweDft* samples;
//    int l;
//    int bgBit;
};

TrgswKey* trgswInitKey(TrlweKey* trlweKey, int l, int bgBit);

void initTrgswSample(Trgsw* trgsw, int l, int bgBit, int k, int N);

void initTrgswDftSample(TrgswDft* trgsw, int l, int bgBit, int k, int N);

void initTrgswDftSample(TrgswDft* trgswDft, int n, int l, int bgBit, int k, int N);

#endif //HLS_YATFHE_TRGSW_H
