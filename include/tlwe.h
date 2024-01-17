//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H

#include "torus.h"
#include "yatfhe_parameters.h"
#include <vector>

struct TlweKey {
//    Integer* bskDft{ new Integer };
    int n {};
    double sigma {};
    std::vector<Integer> s {}; // n
};

void initTlweKey(TlweKey& key, const int n, const double sigma);

void newBinaryTlweKey(TlweKey& key, const YatfheParameters& parameters);

void lweKeyGen(TlweKey& result, int n);

void deleteLweKey(TlweKey& key);

#endif //HLS_YATFHE_TLWE_H
