//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include <vector>
#include "yatfhe_parameters.h"
#include "trlwe.h"

struct Trgsw {
    std::vector<Trlwe> trlweSamples {}; // l * (k + 1)
    std::vector<std::vector<Trlwe>> trlweSamples2 {};
    int l;
//    int bgBit;

    explicit Trgsw(const YatfheParameters& p) :
        trlweSamples(p.l * (p.k + 1), Trlwe(p.k, p.N)),
        trlweSamples2(p.k + 1, std::vector<Trlwe>(p.l, Trlwe(p.k, p.N))),
        l(p.l) {};
};

struct TrgswKey {
    TrlweKey trlweKey {};
    int l {};
    int bgBit {};

    explicit TrgswKey(const YatfheParameters& p) :
        l(p.l),
        bgBit(p.bgBit),
        trlweKey(TrlweKey(p.k, p.N)) {};
};

struct TrgswDft {
//    std::vector<TrlweDft> trlweDftSamples; // l *  (k + 1)
    std::vector<std::vector<TrlweDft>> trlweDftSamples2; // l *  (k + 1)
    int l;
//    int bgBit;

    explicit TrgswDft(const YatfheParameters& p) :
//        trlweDftSamples(p.l * (p.k + 1), TrlweDft(p.k, p.N)),
        trlweDftSamples2(p.k + 1, std::vector<TrlweDft>(p.l, TrlweDft(p.k, p.N))),
        l(p.l) {};
};

//void trgswInitKey(TrgswKey& trgswKey, TrlweKey& trlweKey, const YatfheParameters& param);
//
//void initTrgswSample(Trgsw& trgsw, const YatfheParameters& param);
//
//void initTrgswDftSample(TrgswDft &trgswDftSample, const YatfheParameters& param);

void genNoiselessTrgswSample(Trgsw &trgswSample, Torus msg, const YatfheParameters& param);

//void deleteTrgswKey(TrgswKey& trgswKey);

#endif //HLS_YATFHE_TRGSW_H
