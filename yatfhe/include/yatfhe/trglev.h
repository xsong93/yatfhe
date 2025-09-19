//
// Created by Xintong Song on 2024/4/22.
//

#ifndef HLS_YATFHE_TRGLEV_H
#define HLS_YATFHE_TRGLEV_H

#include "yatfhe/trlwe.h"

struct Trglev {
    std::vector<Trlwe> trlwes; // l
    int l;

    explicit Trglev(const YatfheParameters& param) :
            trlwes(param.l,  Trlwe(param.k, param.N)),
            l(param.l) {};
};

struct TrglevDft {
    std::vector<TrlweDft> trlweDfts; // l
    int l;

    explicit TrglevDft(const YatfheParameters& param) :
            trlweDfts(param.l, TrlweDft(param.k, param.N)),
            l(param.l) {};
};

void encTrglevSingleSample(Trglev& output, const TrlweKey& trlweKey, Torus input, const YatfheParameters& param);

void encTrglevMultiSample(Trglev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void multTrglevWithConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

void multDecomposedTglevWithConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGLEV_H
