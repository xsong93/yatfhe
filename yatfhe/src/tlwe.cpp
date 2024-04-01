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

void lweKeyGen(TlweKey& key, const int n) {
    for (auto i = 0; i < n; i++) {
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
double symDecTlweSample(Tlwe& in, const TlweKey& key, const int torusBase) {
    Torus aXs= 0;
    for (auto i = 0; i < key.n; i++) {
        aXs += in.a[i] * key.s[i];
    }
    return roundError(torus32ToDouble(in.b - aXs), torusBase);
}

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input) {
    const auto newMod = output.mod;
    output.b = modSwitchFromTorus32(input.b, newMod);
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = modSwitchFromTorus32(input.a[i], newMod);
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
