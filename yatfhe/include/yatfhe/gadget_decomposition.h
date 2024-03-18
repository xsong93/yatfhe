//
// Created by Xintong Song on 2024/3/15.
//

#ifndef HLS_YATFHE_GADGET_DECOMPOSITION_H
#define HLS_YATFHE_GADGET_DECOMPOSITION_H

#include <iostream>
#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"

struct DecomposedData {
    std::vector<Integer> value;
    Integer sign {};

    explicit DecomposedData(int size) : sign(1), value(size, 0) {};
};

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

void gadgetDecompose(DecomposedData& out, Integer in, const YatfheParameters& param);

Integer recompose(const DecomposedData& digits, const YatfheParameters& param);

std::vector<Integer> decomposeOverB(Binary in, const YatfheParameters& param);

void signedGadgetDecomposition(DecomposedData& res, Integer in, const YatfheParameters& param);

int genOffset(int radixBits, int bHalf, int l, int torusBits);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_GADGET_DECOMPOSITION_H
