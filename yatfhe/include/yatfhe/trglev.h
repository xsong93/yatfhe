//
// Created by Xintong Song on 2024/4/22.
//

#ifndef HLS_YATFHE_TRGLEV_H
#define HLS_YATFHE_TRGLEV_H

#include "yatfhe/trlwe.h"

struct Trglev {
    std::vector<Trlwe> trlwes; // l
    std::vector<TrlweDft> trlweDfts; // 2l
    int l;
    int lDft;

    explicit Trglev(const YatfheParameters& param) :
            l(param.l),
            lDft(param.l * (param.dftBits / param.torusBits)),
            trlwes(param.l,  Trlwe(param.k, param.N)),
            trlweDfts(param.l * (param.dftBits / param.torusBits), TrlweDft(param.k, param.N)) {};
};

void encTrglevSingleSample(Trglev& output, const TrlweKey& trlweKey, Torus input, const YatfheParameters& param);

void encTrglevMultiSample(Trglev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void multTrglevWithConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

void multDecomposedTglevWithConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGLEV_H
