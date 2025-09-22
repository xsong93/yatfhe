//
// Created by Xintong Song on 2024/4/22.
//

#ifndef HLS_YATFHE_TRGLEV_H
#define HLS_YATFHE_TRGLEV_H

#include "yatfhe/trlwe.h"

struct Trlev {
    std::vector<Trlwe> trlwes; // l
    int l;

    explicit Trlev(const YatfheParameters& param) :
            trlwes(param.l,  Trlwe(param.k, param.N)),
            l(param.l) {}
};

struct TrlevDft {
    std::vector<TrlweDft> trlweDfts; // l
    int l;

    explicit TrlevDft(const YatfheParameters& param) :
            trlweDfts(param.l, TrlweDft(param.k, param.N)),
            l(param.l) {}
};

template<typename TrlevType>
void addTrlev(TrlevType& out, const TrlevType& in1, const TrlevType& in2) {
    const auto L = out.l;
    for (size_t l = 0; l < L; l++) {
        addTrlwe(out.trlwes[l], in1.trlwes[l], in2.trlwes[l]);
    }
}

void encTrlevSingleSample(Trlev& output, const TrlweKey& trlweKey, Torus input, const YatfheParameters& param);

void encTrlevMultiSample(Trlev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void decTrlev(TorusPolynomial& output, const Trlev& input, const TrlweKey& trlweKey, const YatfheParameters& param);

void multTrlevWithConst(Trlwe& output, const Trlev& input, Integer num, const YatfheParameters& param);

void multDecomposedTrlevWithConst(Trlwe& output, const Trlev& input, Integer num, const YatfheParameters& param);

void rotateTrlev(Trlev& trglev, const int rot, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGLEV_H
