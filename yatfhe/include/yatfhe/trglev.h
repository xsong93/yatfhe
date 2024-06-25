//
// Created by Xintong Song on 2024/4/22.
//

#ifndef HLS_YATFHE_TRGLEV_H
#define HLS_YATFHE_TRGLEV_H

#include <iostream>
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

void trglevEncSingleSample(Trglev& output, const TrlweKey& trlweKey, Torus input, const YatfheParameters& param);

void trglevEncMultiSample(Trglev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void trglevMultConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

void trglevDotMultConst(Trlwe& output, const Trglev& input, const std::vector<Integer>& nums, const YatfheParameters& param);

void decomposedTglevMultConst(Trlwe& output, const Trglev& input, Integer num, const YatfheParameters& param);

#endif //HLS_YATFHE_TRGLEV_H
