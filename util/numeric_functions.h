//
// Created by Xintong Song on 2023/12/11.
//

#ifndef HLS_YATFHE_NUMERIC_FUNCTIONS_H
#define HLS_YATFHE_NUMERIC_FUNCTIONS_H

#include <cstdlib>
#include <random>
#include "torus.h"

using namespace std;

extern default_random_engine generator;
static const int64_t twoP32 = INT64_C(1) << 32; // 2^32

Torus gaussian32(Torus message, double sigma);

Torus dToT32(double d);

#endif //HLS_YATFHE_NUMERIC_FUNCTIONS_H
