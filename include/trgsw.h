//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include <vector>
#include "yatfhe_parameters.h"
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

void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const YatfheParameters& yatfheParameters);

void initTrgswSample(Trgsw& trgsw, const YatfheParameters& yatfheParameters);

void initTrgswDftSample(TrgswDft &trgswDftSample, const YatfheParameters& yatfheParameters);

void deleteTrgswKey(TrgswKey& trgswKey);

#endif //HLS_YATFHE_TRGSW_H
