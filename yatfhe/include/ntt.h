//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_H
#define HLS_YATFHE_NTT_H

#include <vector>
#include "torus.h"
#include "polynomial.h"

constexpr uint64_t MODULUS = 0xffffffff00000001UL;

void applyNtt(LagrangePolynomial& lagrangePolynomial, const TorusPolynomial& torusPolynomial);

void applyIntt(LagrangePolynomial& torusPolynomial, LagrangePolynomial& lagrangePolynomial);

void bitRevShuffle(std::vector<uint64_t>& x, int N);

uint64_t modAdd(uint64_t x, uint64_t y);

uint64_t modSub(uint64_t x, uint64_t y);

uint64_t modMul(uint64_t x, uint64_t y);

#endif //HLS_YATFHE_NTT_H
