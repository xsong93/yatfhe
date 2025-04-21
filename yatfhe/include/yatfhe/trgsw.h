//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_TRGSW_H
#define HLS_YATFHE_TRGSW_H

#include <vector>
#include <hexl/hexl.hpp>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trlwe.h"

struct TrgswMP {
    std::vector<std::vector<Trlwe>> c;
    std::vector<Trlwe> cPrime;
    int l;
    int k;

    explicit TrgswMP(const YatfheParameters& p) :
            c(p.l, std::vector<Trlwe>(p.k, Trlwe(p.k, p.N))),
            cPrime(p.l, Trlwe(p.k, p.N)),
            l(p.l),
            k(p.k) {};
};

struct TrgswMPDft {
    std::vector<std::vector<TrlweDft>> c;
    std::vector<TrlweDft> cPrime;
    int l;
    int k;

    explicit TrgswMPDft(const YatfheParameters& p) :
            c(p.l, std::vector<TrlweDft>(p.k, TrlweDft(p.k, p.N))),
            cPrime(p.l, TrlweDft(p.k, p.N)),
            l(p.l),
            k(p.k) {};
};

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

struct Trgsw8 {
    std::vector<std::vector<Trlwe8>> trlweSamples {};
    int l;
    int k;

    explicit Trgsw8(const YatfheParameters& p) :
            trlweSamples(p.l, std::vector<Trlwe8>(p.k + 1, Trlwe8(p.k, p.N))),
            l(p.l),
            k(p.k) {};

    Trgsw8(const int l, const int k, const int N) :
            trlweSamples(l, std::vector<Trlwe8>(k + 1, Trlwe8(k, N))),
            l(l),
            k(k) {};
};

struct TrgswDft {
    std::vector<std::vector<TrlweDft>> trlweDftSamples; // l *  (k + 1)
    int l;
    int k;
//    int bgBit;

    explicit TrgswDft(const YatfheParameters& p) :
//            trlweDftSamples(p.k + 1, std::vector<TrlweDft>(p.l, TrlweDft(p.k, p.N))),
            trlweDftSamples(p.l, std::vector<TrlweDft>(p.k + 1, TrlweDft(p.k, p.N))),
            l(p.l),
            k(p.k) {};
};

struct TrgswDft14 {
    std::vector<std::vector<TrlweDft14>> trlweDftSamples; // l *  (k + 1)
    int l;
    int k;

    explicit TrgswDft14(const YatfheParameters& p) :
            trlweDftSamples(p.l, std::vector<TrlweDft14>(p.k + 1, TrlweDft14(p.k, p.N))),
            l(p.l),
            k(p.k) {};
};

struct TrgswDft24 {
    std::vector<std::vector<TrlweDft24>> trlweDftSamples; // l *  (k + 1)
    int l;
    int k;

    explicit TrgswDft24(const YatfheParameters& p) :
            trlweDftSamples(p.l, std::vector<TrlweDft24>(p.k + 1, TrlweDft24(p.k, p.N))),
            l(p.l),
            k(p.k) {};

    TrgswDft24(const int l, const int k, const int N) :
            trlweDftSamples(l, std::vector<TrlweDft24>(k + 1, TrlweDft24(k, N))),
            l(l),
            k(k) {};
};

struct TrgswKey {
    TrlweKey trlweKey;

    explicit TrgswKey(const YatfheParameters& p) :
        trlweKey(TrlweKey(p)) {};
};

template<typename TrgswType>
void trgswAdd(TrgswType& out, const TrgswType& in1, const TrgswType& in2) {
    auto L = out.l;
    auto K = out.k;
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K + 1; k++) {
            trlweAdd(out.trlweSamples[l][k], in1.trlweSamples[l][k], in2.trlweSamples[l][k]);
        }
    }
}

template<typename TrgswDftType>
void trgswAddNtt(TrgswDftType& out, const TrgswDftType& in1, const TrgswDftType& in2) {
    auto L = out.l;
    auto K = out.k;
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K + 1; k++) {
            trlweAddNtt(out.trlweDftSamples[l][k], in1.trlweDftSamples[l][k], in2.trlweDftSamples[l][k]);
        }
    }
}

template<typename TrgswType>
void trgswSub(TrgswType& out, const TrgswType& in1, const TrgswType& in2) {
    auto L = out.l;
    auto K = out.k;
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K + 1; k++) {
            trlweSub(out.trlweSamples[l][k], in1.trlweSamples[l][k], in2.trlweSamples[l][k]);
        }
    }
}

template<typename TrgswDftType>
void trgswSubNtt(TrgswDftType& out, const TrgswDftType& in1, const TrgswDftType& in2) {
    auto L = out.l;
    auto K = out.k;
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K + 1; k++) {
            trlweSubNtt(out.trlweDftSamples[l][k], in1.trlweDftSamples[l][k], in2.trlweDftSamples[l][k]);
        }
    }
}

void trgswRotate(Trgsw& trgsw, int rot, const YatfheParameters& param);

void trgswRotateNtt(TrgswDft& trgswDft, const int rot, const YatfheParameters& param);

void trgswMPEncrypt(TrgswMP& trgswMP, Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param);

void trgswMPEncryptNtt(TrgswMP& trgswMP, TrgswMPDft& trgswMPDft, Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param);

void trgswMPEncryptLow(TrgswMP& trgswMP, Integer mu, const TrgswKey& trgswKey, const YatfheParameters& param);

void trgswMPEncryptLowNtt(TrgswMP& trgswMP, TrgswMPDft& trgswMPDft, Integer mu, const TrgswKey& trgswKey, const YatfheParameters& param);

void trgswEncZero(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswAddInteger(Trgsw& trgsw, Integer mu, const int pos, const YatfheParameters& param);

void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, Integer mu, const int pos, const YatfheParameters& param);

void trgswEncrypt(Trgsw& trgsw, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param);

void trgswEncryptNtt(Trgsw& trgsw, TrgswDft& trgswDft, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param);

void trgswEncryptApproxCRT(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

Integer trgswDecrypt(const Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey);

Integer trgswDecryptNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void trgswMCRTDecomp(std::vector<Trgsw8>& out, const Trgsw& in, const YatfheParameters& param);

void trgswMCRTToCRT(std::vector<Trgsw8>& trgsw, const YatfheParameters& param);

void trgswCRTRecomp(Trgsw& out, const std::vector<Trgsw8>& in, const YatfheParameters& param);

void trgswExternalProduct(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswExternalProductApproxCRT(std::vector<Trlwe8>& output, const std::vector<Trgsw8>& trgswInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param);

void trgswExternalProductApproxCRTNtt(std::vector<Trlwe8>& output, const std::vector<TrgswDft24>& trgswDftInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param);

void trgswMPExternalProduct(Trlwe& output, const TrgswMP& trgswMPInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswMPExternalProductNtt(Trlwe& output, const TrgswMPDft& trgswMPInput, const Trlwe& trlweInput, const YatfheParameters& param);

void trgswMPExternalProductDecomp(DecomposedTrlwe& output, const TrgswMP& trgswMPInput, const DecomposedTrlwe& trlweInput, const YatfheParameters& param);

void trgswMPExternalProductDecompNtt(DecomposedTrlweDft& output, const TrgswMPDft& trgswMPInput, const DecomposedTrlweDft& trlweInput, const YatfheParameters& param);

void trgswMPInternalProduct(TrgswMP& output, const TrgswMP& input1, const TrgswMP& input2, const YatfheParameters& param);

void trgswMPInternalProductNtt(TrgswMP& output, const TrgswMP& input1, const TrgswMPDft& input2, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGSW_H
