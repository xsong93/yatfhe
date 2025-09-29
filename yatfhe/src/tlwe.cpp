//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yautil/control_helper.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"

using namespace std;

void genTlweKey(TlweKey& key) {
    for (auto i = 0; i < key.n; i++) {
#ifdef TERNRY
        key.s[i] = ternaryDistrib(rng);
#else
        key.s[i] = binaryDistrib(rng);
#endif
    }
}

// b = aj * sj + u + e
void symEncTlwe(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorusDistrib()(rng);
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(tlweSample.a[i]) * key.s[i];
        }
    }
    Torus muE = addGaussianNoise(message, key.sigma);
    tlweSample.b = addTorus(muE, static_cast<Torus>(longModP(tmp, TORUS_Q)));
}

// mu = b - as
double symDecTlweToDouble(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, TORUS_Q));
    return roundError(torus32ToDouble(subTorus(in.b, aXs)), torusBase);
}

Torus symDecTlweToTorus(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, TORUS_Q));
    return roundTorusError(subTorus(in.b, aXs), torusBase);
}

Integer symDecTlweToInt(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, TORUS_Q));
    return modSwitchFromTorus32(roundTorusError(subTorus(in.b, aXs), torusBase), torusBase);
}

Torus calTlweError(Tlwe& in, const TlweKey& key, Torus mu) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, TORUS_Q));
    return subTorus(in.b, aXs) - mu;
}

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input) {
    const auto newMod = output.mod;
    output.b = static_cast<int32_t>(modSwitchFromTorus32(input.b, newMod));
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = static_cast<int32_t>(modSwitchFromTorus32(input.a[i], newMod));
    }
}

void addTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = addTorus(input1.a[i], input2.a[i]);
    }
    output.b = addTorus(input1.b, input2.b);
}

void subTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = subTorus(input1.a[i], input2.a[i]);
    }
    output.b = subTorus(input1.b, input2.b);
}

// output -= input
void subTlweInPlace(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = subTorus(output.a[i], input.a[i]);
    }
    output.b = subTorus(output.b, input.b);
}

void copyTlwe(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = input.a[i];
    }
    output.b = input.b;
}

void resetTlweToZero(Tlwe& tlwe) {
    tlwe.b = 0;
    for (auto& ai : tlwe.a) {
        ai = 0;
    }
}
