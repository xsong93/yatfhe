//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

void lweKeyGen(TlweKey& key) {
    for (auto i = 0; i < key.n; i++) {
        key.s[i] = binaryDistrib(rng);
    }
}

// b = aj * sj + u + e
void symEncTlweSample(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorusDistrib(rng);
        if (key.s[i] != 0) {
            tmp += + key.s[i] * static_cast<int64_t>(tlweSample.a[i]);
        }
    }
    tlweSample.b = addGaussianNoise(message, key.sigma) + static_cast<Torus>(longModP(tmp, TORUS_Q));
}

// mu = b - as
double symDecTlweSampleToDouble(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = longModP(tmp, TORUS_Q);
    int64_t tmp2 = static_cast<int64_t>(in.b) - aXs;
    auto bMinusAxS = static_cast<Torus>(longModP(tmp2, TORUS_Q));
    return roundError(torus32ToDouble(bMinusAxS), torusBase);
}

Torus symDecTlweSampleToTorus(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto tmp2 = static_cast<int64_t>(in.b) - longModP(tmp, TORUS_Q);
    auto bMinusAxS = static_cast<Torus>(longModP(tmp2, TORUS_Q));
    return roundTorusError(bMinusAxS, torusBase);
}

Integer symDecTlweSampleToInt(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto tmp2 = static_cast<int64_t>(in.b) - longModP(tmp, TORUS_Q);
    auto bMinusAxS = static_cast<Torus>(longModP(tmp2, TORUS_Q));
    return modSwitchFromTorus32(roundTorusError(bMinusAxS, torusBase), torusBase);
}

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input) {
    const auto newMod = output.mod;
    output.b = static_cast<int32_t>(modSwitchFromTorus32Pos(input.b, newMod));
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = static_cast<int32_t>(modSwitchFromTorus32Pos(input.a[i], newMod));
    }
}

void lweAdd(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    int64_t tmp = 0;
    for (auto i = 0; i < output.n; i++) {
        tmp = static_cast<int64_t>(input1.a[i]) + static_cast<int64_t>(input2.a[i]);
        output.a[i] = static_cast<Torus>(longModP(tmp, TORUS_Q));
    }
    tmp = static_cast<int64_t>(input1.b) + static_cast<int64_t>(input2.b);
    output.b = static_cast<Torus>(longModP(tmp, TORUS_Q));
}

void lweSub(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    int64_t tmp = 0;
    for (auto i = 0; i < output.n; i++) {
        tmp = static_cast<int64_t>(input1.a[i]) - static_cast<int64_t>(input2.a[i]);
        output.a[i] = static_cast<Torus>(longModP(tmp, TORUS_Q));
    }
    tmp = static_cast<int64_t>(input1.b) - static_cast<int64_t>(input2.b);
    output.b = static_cast<Torus>(longModP(tmp, TORUS_Q));
}

// output -= input
void lweSubTo(Tlwe& output, const Tlwe& input) {
    int64_t tmp = 0;
    for (auto i = 0; i < output.n; i++) {
        tmp = static_cast<int64_t>(output.a[i]) - static_cast<int64_t>(input.a[i]);
        output.a[i] = static_cast<Torus>(longModP(tmp, TORUS_Q));
    }
    tmp = static_cast<int64_t>(output.b) - static_cast<int64_t>(input.b);
    output.b = static_cast<Torus>(longModP(tmp, TORUS_Q));
}

void tlweCopy(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = input.a[i];
    }
    output.b = input.b;
}

void tlweCLear(Tlwe& tlwe) {
    tlwe.b = 0;
    for (auto& ai : tlwe.a) {
        ai = 0;
    }
}
