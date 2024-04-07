//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include <vector>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trlwe.h"

struct Trgsw {
    std::vector<std::vector<Trlwe>> trlweSamples {};
    int l;
//    int bgBit;

    explicit Trgsw(const YatfheParameters& p) :
//            trlweSamples(p.k + 1, std::vector<Trlwe>(p.l, Trlwe(p.k, p.N))),
            trlweSamples(p.lDft, std::vector<Trlwe>(p.k + 1, Trlwe(p.k, p.N))),
            l(p.lDft) {};
};

struct TrgswDft {
    std::vector<std::vector<TrlweDft>> trlweDftSamples; // l *  (k + 1)
    int l;
//    int bgBit;

    explicit TrgswDft(const YatfheParameters& p) :
//            trlweDftSamples(p.k + 1, std::vector<TrlweDft>(p.l, TrlweDft(p.k, p.N))),
            trlweDftSamples(p.lDft, std::vector<TrlweDft>(p.k + 1, TrlweDft(p.k, p.N))),
            l(p.lDft) {};
};

struct TrgswKey {
    TrlweKey trlweKey;

    explicit TrgswKey(const YatfheParameters& p) :
        trlweKey(TrlweKey(p.k, p.N)) {};
};

void trgswEncrypt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey, Integer mu);

Integer trgswDecrypt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey);

void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, int64_t mu, const YatfheParameters& param);

void trgswExternalProduct(Trlwe& output, TrgswDft& trgswInput, Trlwe& trlweInput, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGSW_H
