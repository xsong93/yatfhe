//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include "trlwe.h"

struct Trgsw {
    Trlwe* trlweSamples{new Trlwe }; // l * (k + 1)
//    int l;
//    int bgBit;
};

struct TrgswKey {
    TrlweKey* trlweKey{ new TrlweKey };
    int l{};
    int bgBit{};
};

struct TrgswDft {
    TrlweDft* trlweDftSamples{new TrlweDft };
//    int l;
//    int bgBit;
};

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, int l, int bgBit);

void initTrgswSample(Trgsw& trgsw, int l, int bgBit, int k, int N);

void initTrgswDftSample(TrgswDft& trgswDft, int l, int bgBit, int k, int N);

void deleteTrgswKey(TrgswKey& trgswKey);

#endif //HLS_YATFHE_TRGSW_H
