//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H


//#include <stdio.h>
//#include <stdint.h>
//#include <string.h>
//#include <math.h>
//#include <assert.h>
//#include <string.h>
//#include <stdbool.h>
#include "torus.h"

struct TlweKey {
    Integer * s;
    int n;
    double sigma;
};

void tlweInitKey(TlweKey *key, int n, double sigma);

//TlweKey* tlweInitKey(int n, double sigma);

TlweKey *tlweNewBinaryKey(int n, double sigma);

void lweKeyGen(TlweKey* result, int n);

void deleteLweKey(TlweKey* key);

#endif //HLS_YATFHE_TLWE_H
