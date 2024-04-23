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
            trlweSamples(p.l, std::vector<Trlwe>(p.k + 1, Trlwe(p.k, p.N))),
            l(p.l) {};
};

struct DecomposedTrgsw {
    std::vector<Trgsw> trgswSamples {};
    int l;

    explicit DecomposedTrgsw(const YatfheParameters& p) :
//            trlweSamples(p.k + 1, std::vector<Trlwe>(p.l, Trlwe(p.k, p.N))),
            trgswSamples(p.l2, Trgsw(p)),
            l(p.l2) {};
};

struct TrgswDft {
    std::vector<std::vector<TrlweDft>> trlweDftSamples; // l *  (k + 1)
    int l;
//    int bgBit;

    explicit TrgswDft(const YatfheParameters& p) :
//            trlweDftSamples(p.k + 1, std::vector<TrlweDft>(p.l, TrlweDft(p.k, p.N))),
            trlweDftSamples(p.l, std::vector<TrlweDft>(p.k + 1, TrlweDft(p.k, p.N))),
            l(p.l) {};
};

struct TrgswKey {
    TrlweKey trlweKey;

    explicit TrgswKey(const YatfheParameters& p) :
        trlweKey(TrlweKey(p)) {};
};

void trgswEncrypt(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

void trgswEncryptNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

Integer trgswDecrypt(const Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey);

Integer trgswDecryptNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswExternalProduct(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGSW_H
