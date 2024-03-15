//
// Created by Xintong Song on 2024/3/15.
//

#ifndef HLS_YATFHE_GADGET_DECOMPOSITION_H
#define HLS_YATFHE_GADGET_DECOMPOSITION_H

#include <iostream>
#include "yatfhe/torus.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"

std::vector<Torus> genGadgetVector(int radixBits, int l, int torusBits);

UnsignedInteger recompose(const std::vector<UnsignedInteger>& digits, const YatfheParameters& param);

std::vector<Integer> decomposeOverB(const Integer in, const YatfheParameters& param);

void signedGadgetDecomposition(std::vector<Torus>& res, const Torus input, const YatfheParameters& param);

int genOffset(int radixBits, int bHalf, int l, int torusBits);

void gadgetDecomposeTrlwe(DecomposedTrlwe& output, Trlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_GADGET_DECOMPOSITION_H
