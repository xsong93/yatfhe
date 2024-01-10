//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H

#include "torus.h"
#include <vector>

struct TlweKey {
//    Integer* bskDft{ new Integer };
    int n {};
    double sigma {};
    std::vector<Integer> s {}; // n
};

void tlweInitKey(TlweKey& key, int n, double sigma);

void tlweNewBinaryKey(TlweKey& key, int n, double sigma);

void lweKeyGen(TlweKey& result, int n);

void deleteLweKey(TlweKey& key);

#endif //HLS_YATFHE_TLWE_H
