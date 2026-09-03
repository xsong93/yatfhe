//
// Created by Xintong Song on 2024/4/22.
//

#ifndef YATFHE_TRGLEV_H
#define YATFHE_TRGLEV_H

#include "yatfhe/trlwe.h"

struct Trlev {
    std::vector<Trlwe> trlwes; // l
    int l{};

    Trlev() = default;

    explicit Trlev(const YatfheParameters& param) :
            trlwes(param.l,  Trlwe(param.k, param.N)),
            l(param.l) {}
    Trlev(const YatfheParameters& param, const int level) :
            trlwes(level,  Trlwe(param.k, param.N)),
            l(level) {}
};

struct TrlevDft {
    std::vector<TrlweDft> trlweDfts; // l
    int l{};

    TrlevDft() = default;

    explicit TrlevDft(const YatfheParameters& param) :
            trlweDfts(param.l, TrlweDft(param.k, param.N)),
            l(param.l) {}
    TrlevDft(const YatfheParameters& param, const int level) :
            trlweDfts(level, TrlweDft(param.k, param.N)),
            l(level) {}
};

template<typename TrlevType>
void addTrlev(TrlevType& out, const TrlevType& in1, const TrlevType& in2) {
    const auto L = out.l;
    for (size_t l = 0; l < L; l++) {
        addTrlwe(out.trlwes[l], in1.trlwes[l], in2.trlwes[l]);
    }
}

template<typename TrlevType>
void subTrlev(TrlevType& out, const TrlevType& in1, const TrlevType& in2) {
    const auto L = out.l;
    for (size_t l = 0; l < L; l++) {
        subTrlwe(out.trlwes[l], in1.trlwes[l], in2.trlwes[l]);
    }
}

void encTrlevSingleSample(Trlev& output, const TrlweKey& trlweKey, Torus input, int pos, const YatfheParameters& param);

// Encode a scalar `input` as monomial X^{monomialIndex} in the RLWE plaintext polynomial domain.
void encTrlevSingleSampleMonomial(Trlev& output, const TrlweKey& trlweKey, Torus input, int monomialIndex,
                                   const YatfheParameters& param);

void encTrlevMultiSample(Trlev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void symEncTrlevWithKey(Trlev& output, const TrlweKey& trlweKey, const vector<TorusPolynomial>& inputs, bool isPos, const YatfheParameters& param);

void symEncTrlevWithKeyNtt(TrlevDft& output, const TrlweKey& trlweKey, const vector<IntPolynomial>& inputs, bool isPos, const YatfheParameters& param);

void decTrlev(TorusPolynomial& output, const Trlev& input, const TrlweKey& trlweKey, const YatfheParameters& param);

void multTrlevWithConst(Trlwe& output, const Trlev& input, Torus num, const YatfheParameters& param);

void multDecomposedTrlevWithConst(Trlwe& output, const Trlev& input, Torus num, const YatfheParameters& param);

void rotateTrlev(Trlev& trglev, const int rot, const YatfheParameters& param);

#endif //YATFHE_TRGLEV_H
