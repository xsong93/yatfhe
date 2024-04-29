//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"

using namespace std;

void lweKeyGen(TlweKey& key) {
    for (auto i = 0; i < key.n; i++) {
        key.s[i] = binaryDistrib(rng);
    }
//    printArray(key.s, "TlweKey");
}

// b = aj * sj + u + e
void symEncTlweSample(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    tlweSample.b = addGaussianNoise(message, key.sigma); // error term
    for (auto i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorusDistrib(rng);
        tlweSample.b += key.s[i] * tlweSample.a[i];
    }
}

// mu = b - as
double symDecTlweSampleToDouble(Tlwe& in, const TlweKey& key, const int torusBase) {
    Torus aXs = 0;
    for (auto i = 0; i < key.n; i++) {
        aXs += in.a[i] * key.s[i];
    }
    return roundError(torus32ToDouble(in.b - aXs), torusBase);
}

Torus symDecTlweSampleToTorus(Tlwe& in, const TlweKey& key, const int torusBase) {
    Torus aXs = 0;
    for (auto i = 0; i < key.n; i++) {
        aXs += in.a[i] * key.s[i];
    }
    return roundTorusError(in.b - aXs, torusBase);
}

Integer symDecTlweSampleToInt(Tlwe& in, const TlweKey& key, const int torusBase) {
    Torus aXs = 0;
    for (auto i = 0; i < key.n; i++) {
        aXs += in.a[i] * key.s[i];
    }
    return modSwitchFromTorus32(roundTorusError(in.b - aXs, torusBase), torusBase);
}

//todo: nt
void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input) {
    const auto newMod = output.mod;
    output.b = modSwitchFromTorus32Pos(input.b, newMod);
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = modSwitchFromTorus32Pos(input.a[i], newMod);
    }
}

void lweAdd(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = input1.a[i] + input2.a[i];
    }
    output.b = input1.b + input2.b;
}

void lweSub(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = input1.a[i] - input2.a[i];
    }
    output.b = input1.b - input2.b;
}

// output -= input
void lweSubTo(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] -= input.a[i];
    }
    output.b -= input.b;
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
