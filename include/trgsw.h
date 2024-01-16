//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include <vector>
#include "trlwe.h"

struct Trgsw {
//    Trlwe* trlweSamples{new Trlwe }; // l * (k + 1)
    std::vector<Trlwe> trlweSamples {}; // l * (k + 1)
//    int l;
//    int bgBit;
};

struct TrgswKey {
//    TrlweKey* trlweKey{ new TrlweKey };
    TrlweKey trlweKey {};
    int l {};
    int bgBit {};
};

struct TrgswDft {
//    TrlweDft* trlweDftSamples{new TrlweDft };
    std::vector<TrlweDft> trlweDftSamples {};
//    int l;
//    int bgBit;
};

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, int l, int bgBit);

void initTrgswSample(Trgsw& trgsw, const TrgswKey& trgswKey);

void initTrgswDftSample(TrgswDft& trgswDft, const TrgswKey& trgswKey);

void allocNewTrgswDftSample(TrgswDft &trgswDftSample, int l, int Bg_bit, int k, int N);

void deleteTrgswKey(TrgswKey& trgswKey);

#endif //HLS_YATFHE_TRGSW_H
