//
// Created by Xintong Song on 2024/3/15.
//

#ifndef HLS_YATFHE_GADGET_DECOMPOSITION_H
#define HLS_YATFHE_GADGET_DECOMPOSITION_H

#include <iostream>
#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"

struct DecomposedData {
    std::vector<Integer> value; // l
    int l {};
    Integer sign {1}; // set default to 1 as positive sign

    explicit DecomposedData(int size) : l(size), value(size, 0) {};
};

struct DecomposedDataDft {
    std::vector<NttType> value; // l
    int l {};

    explicit DecomposedDataDft(int size) : l(size), value(size, 0) {};
};

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

void gadgetDecompose(DecomposedData& out, Integer in, const YatfheParameters& param);

Integer recomposeSelf(const DecomposedData& digits, const YatfheParameters& param);

void recomposeFirstHalf(DecomposedData& output, const DecomposedData& lhs, const std::vector<DecomposedData>& mid);

Integer recomposeTwoParts(const DecomposedData& lhs, const std::vector<Integer>& rhs);

void decomposeOverB(std::vector<Integer>& output, Integer in, const YatfheParameters& param);

void signedGadgetDecomposition(DecomposedData& res, Integer in, const YatfheParameters& param);

void signedGadgetDecompositionNtt(DecomposedDataDft& out, NttType in, const YatfheParameters& param);

int genOffset(int radixBits, int bHalf, int l, int torusBits);

#endif //HLS_YATFHE_GADGET_DECOMPOSITION_H
