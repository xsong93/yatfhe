//
// Created for anonymous review.
//

#ifndef HLS_YATFHE_TLWE_H
#define HLS_YATFHE_TLWE_H

#include <vector>

#include "yatfhe_parameters.h"
#include "yatfhe/torus.h"

struct Tlwe {
    std::vector<Torus> a {}; // n
    Torus b {};
    int n {};
    size_t bytes {};

    explicit Tlwe(int n) : a(n), b(0), n(n), bytes(sizeof(Torus) * (n+1)) {};
};

struct ScaledTlwe {
    std::vector<int32_t> a {}; // n
    int32_t b {};
    int n {};
    int mod {};

    ScaledTlwe() = default;

    ScaledTlwe(int mod, int n) : a(n), b(0), n(n), mod(mod) {};
};

struct TlweKey {
    int n {};
    double sigma {};
    std::vector<Binary> s {}; // n

    TlweKey(int n, double sigma) : n(n), sigma(sigma), s(n) {};

    explicit TlweKey(int n) : n(n), sigma(0), s(n) {};

    explicit TlweKey(const YatfheParameters& p) : n(p.n), sigma(p.lweStdDev), s(p.n) {};
};

void genTlweKey(TlweKey& key);

void rescaleTlweFromTorus(Tlwe& output, const Tlwe& input);

void rescaleTlweToNewMod(ScaledTlwe& output, const Tlwe& input);

void symEncTlwe(Tlwe& tlweSample, Torus message, const TlweKey& key);

Torus symDecTlweToTorus(Tlwe& in, const TlweKey& key, int torusBase);

Integer symDecTlweToInt(Tlwe& in, const TlweKey& key, int torusBase);

Torus calTlweError(Tlwe& in, const TlweKey& key, Torus mu);

void addTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void subTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void subTlweInPlace(Tlwe& output, const Tlwe& input);

void copyTlwe(Tlwe& output, const Tlwe& input);

void resetTlweToZero(Tlwe& tlwe);

#endif //HLS_YATFHE_TLWE_H
