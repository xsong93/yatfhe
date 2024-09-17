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
    int k;
//    int bgBit;

    explicit Trgsw(const YatfheParameters& p) :
//            trlweSamples(p.k + 1, std::vector<Trlwe>(p.l, Trlwe(p.k, p.N))),
            trlweSamples(p.l, std::vector<Trlwe>(p.k + 1, Trlwe(p.k, p.N))),
            l(p.l),
            k(p.k) {};
};

struct Trgsw16 {
    std::vector<std::vector<Trlwe16>> trlweSamples {};
    int l;
    int k;

    explicit Trgsw16(const YatfheParameters& p) :
            trlweSamples(p.l, std::vector<Trlwe16>(p.k + 1, Trlwe16(p.k, p.N))),
            l(p.l),
            k(p.k) {};
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

struct TrgswDft14 {
    std::vector<std::vector<TrlweDft14>> trlweDftSamples; // l *  (k + 1)
    int l;

    explicit TrgswDft14(const YatfheParameters& p) :
            trlweDftSamples(p.l, std::vector<TrlweDft14>(p.k + 1, TrlweDft14(p.k, p.N))),
            l(p.l) {};
};

struct TrgswDft24 {
    std::vector<std::vector<TrlweDft24>> trlweDftSamples; // l *  (k + 1)
    int l;

    explicit TrgswDft24(const YatfheParameters& p) :
            trlweDftSamples(p.l, std::vector<TrlweDft24>(p.k + 1, TrlweDft24(p.k, p.N))),
            l(p.l) {};
};

struct TrgswKey {
    TrlweKey trlweKey;

    explicit TrgswKey(const YatfheParameters& p) :
        trlweKey(TrlweKey(p)) {};
};

void trgswEncZero(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswAddInteger(Trgsw& trgsw, Integer mu, const YatfheParameters& param);

void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, Integer mu, const YatfheParameters& param);

void trgswEncrypt(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

void trgswEncryptNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

Integer trgswDecrypt(const Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey);

Integer trgswDecryptNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswExternalProduct(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswInput, Trlwe& trlweInput, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGSW_H
