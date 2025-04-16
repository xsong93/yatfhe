//
// Created by Xintong on 25-4-16.
//

#ifndef HLS_YATFHE_NTT_HEXL_H
#define HLS_YATFHE_NTT_HEXL_H

#include <hexl/hexl.hpp>
#include "yatfhe/polynomial.h"


using namespace intel::hexl;

NTT& nttHexl();

void initNttHexl(uint64_t degree, uint64_t q);

void printHexlParams();

void applyNttHexl(LagrangePolynomial& out, const TorusPolynomial& in);

void applyInttHexl(TorusPolynomial& out, const LagrangePolynomial& in);

#endif //HLS_YATFHE_NTT_HEXL_H
