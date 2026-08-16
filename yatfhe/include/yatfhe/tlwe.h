//
// Created by Xintong Song on 2023/12/8.
//

#ifndef YATFHE_TLWE_H
#define YATFHE_TLWE_H

#include <vector>

#include "yatfhe_parameters.h"
#include "yatfhe/torus.h"

struct Tlwe {
    std::vector<LweTorus> a {}; // n
    LweTorus b {};
    int n {};
    size_t bytes {};

    Tlwe() = default;

    explicit Tlwe(int n) : a(n), b(0), n(n), bytes(sizeof(LweTorus) * (n+1)) {};
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
    int errorB {};
    std::vector<Binary> s {}; // n

    TlweKey(int n, int errorB) : n(n), errorB(errorB), s(n) {};

    explicit TlweKey(int n) : n(n), errorB(0), s(n) {};

    explicit TlweKey(const YatfheParameters& p) : n(p.n), errorB(p.lweNoiseB), s(p.n) {};
};

void genTlweKey(TlweKey& key);

void rescaleTlweFromTorus(Tlwe& output, const Tlwe& input);

void rescaleTlweToNewMod(ScaledTlwe& output, const Tlwe& input);

void symEncTlwe(Tlwe& tlweSample, Torus message, const TlweKey& key);

Torus symDecTlweToTorus(Tlwe& in, const TlweKey& key, int torusBase);

Integer symDecTlweToInt(Tlwe& in, const TlweKey& key, int torusBase);

Torus calTlweError(Tlwe& in, const TlweKey& key, Integer mu);

void addTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void subTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2);

void subTlweInPlace(Tlwe& output, const Tlwe& input);

void multTlwe(Tlwe& output, int scalar);

void copyTlwe(Tlwe& output, const Tlwe& input);

void resetTlweToZero(Tlwe& tlwe);

void inverseTlwe(Tlwe& tlwe);

#endif //YATFHE_TLWE_H
