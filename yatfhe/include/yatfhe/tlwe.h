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
    size_t bytes {};

    explicit Tlwe(int n) : n(n), a(n), b(0), bytes(sizeof(Torus) * (n+1)) {};
};

struct ScaledTlwe {
    std::vector<int32_t> a {}; // n
    int32_t b {};
    int n {};
    int mod {};

    ScaledTlwe(int mod, int n) : mod(mod), n(n), a(n), b(0) {};
};

struct TlweKey {
    int n {};
    double sigma {};
    std::vector<Binary> s {}; // n

    TlweKey(int n, double sigma) : n(n), s(n), sigma(sigma) {};

    explicit TlweKey(int n) : n(n), s(n), sigma(0) {};
};

void lweKeyGen(TlweKey& key);

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input);

void symEncTlweSample(Tlwe& tlweSample, Torus message, const TlweKey& key);

double symDecTlweSampleToDouble(Tlwe& in, const TlweKey& key, int torusBase);

Torus symDecTlweSampleToTorus(Tlwe& in, const TlweKey& key, int torusBase);

Integer symDecTlweSampleToInt(Tlwe& in, const TlweKey& key, int torusBase);

void lweAdd(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void lweSub(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void lweSubTo(Tlwe& output, const Tlwe& input);

void tlweCopy(Tlwe& output, const Tlwe& input);

void tlweCLear(Tlwe& tlwe);

#endif //HLS_YATFHE_TLWE_H
