//
// Created by Xintong Song on 2024/6/12.
//

#ifndef HLS_YATFHE_TGLEV_H
#define HLS_YATFHE_TGLEV_H

#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"

struct Tlev {
    std::vector<Tlwe> tlwes; // l
    int l;

    explicit Tlev(const int l, const int n) :
            l(l),
            tlwes(l,  Tlwe(n)) {};
};

void encTglev(Tlev& output, const TlweKey& tlweKey, Torus input, const YatfheParameters& param);

void multTglevWithConst(Tlwe& output, const Tlev& input, Integer num, const YatfheParameters& param);

// void multDecomposedTglevWithConst(Tlwe& output, const Tglev& input, Integer num, const YatfheParameters& param);

#endif //HLS_YATFHE_TGLEV_H
