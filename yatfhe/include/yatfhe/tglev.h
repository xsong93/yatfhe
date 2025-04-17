//
// Created by Xintong Song on 2024/6/12.
//

#ifndef HLS_YATFHE_TGLEV_H
#define HLS_YATFHE_TGLEV_H

#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"

struct Tglev {
    std::vector<Tlwe> tlwes; // l
    int l;

    explicit Tglev(const int l, const int n) :
            l(l),
            tlwes(l,  Tlwe(n)) {};
};

void tglevEnc(Tglev& output, const TlweKey& tlweKey, Torus input, const YatfheParameters& param);

void tglevMultConst(Tlwe& output, const Tglev& input, Integer num, const YatfheParameters& param);

void decomposedTlevMultConst(Tlwe& output, const Tglev& input, Integer num, const YatfheParameters& param);

#endif //HLS_YATFHE_TGLEV_H
