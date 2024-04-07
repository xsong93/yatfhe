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
    std::vector<Integer> value; // l
    int l {};
    Integer sign {};

    explicit DecomposedData(int size) : sign(1), l(size), value(size, 0) {};
};

struct DecomposedDataDft {
    std::vector<uint64_t> value; // l
    int l {};

    explicit DecomposedDataDft(int size) : l(size), value(size, 0) {};
};

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

void gadgetDecompose(DecomposedData& out, Integer in, const YatfheParameters& param);

void gadgetDecomposeNtt(DecomposedDataDft& out, uint64_t in, const YatfheParameters& param);

Integer recompose(const DecomposedData& digits, const YatfheParameters& param);

std::vector<Integer> decomposeOverB(Integer in, const YatfheParameters& param);

void signedGadgetDecomposition(DecomposedData& res, Integer in, const YatfheParameters& param);

void signedGadgetDecompositionNtt(DecomposedDataDft& out, uint64_t in, const YatfheParameters& param);

int genOffset(int radixBits, int bHalf, int l, int torusBits);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param);

void gadgetDecomposeTrlweNtt(DecomposedTrlwe& output, TrlweDft& input, const YatfheParameters& param);

void recomposeTrlwe(Trlwe& output, DecomposedTrlwe& input, const YatfheParameters& param);

void recomposeTrlweNtt(TrlweDft& output, DecomposedTrlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_GADGET_DECOMPOSITION_H
