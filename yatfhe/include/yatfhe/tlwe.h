//
// Created by Xintong Song on 2023/12/8.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H

#include <vector>
#include "yatfhe/torus.h"
#include "yatfhe/yatfhe_parameters.h"

struct Tlwe {
    std::vector<Torus> a {}; // n
    Torus b {};
    int n {};

    explicit Tlwe(int n) : n(n), a(n), b(0) {};
};

struct ScaledTlwe {
    std::vector<int32_t> a {}; // n
    int32_t b {};
    int n {};
    int mod {};

    explicit ScaledTlwe(int mod, int n) : mod(mod), n(n), a(n), b(0) {};
};

struct TlweKey {
    int n {};
    double sigma {};
    std::vector<Integer> s {}; // n

    explicit TlweKey(int n, double sigma) : n(n), s(n), sigma(sigma) {};

    explicit TlweKey(int n) : n(n), s(n), sigma(0) {};
};

struct TlweKeySwitchingKey {
    //todo
};

void lweKeyGen(TlweKey& result, int n);

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input);

void symEncTlweSample(Tlwe& tlweSample, Torus message, const TlweKey& key);

double symDecTlweSample(Tlwe& in, TlweKey& key);

void lweKeySwitch(Tlwe& output, Tlwe& keySwitchingKey, Tlwe& input, YatfheParameters& param);

void lweSubTo(Tlwe& output, Tlwe& input);

void tlweCopy(Tlwe& output, Tlwe& input);

#endif //HLS_YATFHE_TLWE_H
