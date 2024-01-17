//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H

#include "torus.h"
#include "yatfhe_parameters.h"
#include <vector>

struct Tlwe {
    std::vector<Torus> a {}; // n
    Torus b {};
    int n {};

    explicit Tlwe(int n) : n(n), a(n), b(0) {};
};

struct TlweKey {
    int n {};
    double sigma {};
    std::vector<Integer> s {}; // n

    explicit TlweKey(int n, double sigma) : n(n), s(n), sigma(sigma) {};

    explicit TlweKey(int n) : n(n), s(n), sigma(0) {};
};

void initTlweKey(TlweKey& key, const int n, const double sigma);

void initTlweSample(Tlwe& tlwe, int n);

void newBinaryTlweKey(TlweKey& key, const YatfheParameters& parameters);

void lweKeyGen(TlweKey& result, int n);

void symEncryptToTlweSample(Tlwe& tlweSample, Torus message, const TlweKey& key);

void deleteLweKey(TlweKey& key);

#endif //HLS_YATFHE_TLWE_H
