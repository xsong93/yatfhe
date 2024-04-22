//
// Created by Xintong Song on 2024/4/22.
//

#ifndef HLS_YATFHE_TGLEV_H
#define HLS_YATFHE_TGLEV_H

#include <iostream>
#include "yatfhe/trlwe.h"

struct Tglev {
    std::vector<Trlwe> trlwes; // l
    std::vector<TrlweDft> trlweDfts; // 2l
    int l;
    int lDft;

    explicit Tglev(YatfheParameters param) :
            l(param.l),
            lDft(param.l * (param.dftBits / param.torusBits)),
            trlwes(param.l,  Trlwe(param.k, param.N)),
            trlweDfts(param.l * (param.dftBits / param.torusBits), TrlweDft(param.k, param.N)) {};
};

void tglevEncSingleSample(Tglev& output, const TrlweKey& trlweKey, Torus input, const YatfheParameters& param);

void tglevEncMultiSample(Tglev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param);

void tglevMultConst(Trlwe& output, const Tglev& input, int num, const YatfheParameters& param);

void decomposedTglevMultConst(Trlwe& output, const Tglev& input, int num, const YatfheParameters& param);

#endif //HLS_YATFHE_TGLEV_H
