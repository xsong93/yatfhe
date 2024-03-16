//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

void lweKeyGen(TlweKey& key, const int n) {
    uniform_int_distribution<Binary> distribution(0, 1);
    for (int i = 0; i < n; i++) {
        key.s[i] = distribution(rng);
    }
//    printArray(key.s, "TlweKey");
}

// b = aj * sj + u + e
void symEncTlweSample(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    tlweSample.b = addGaussianNoise(message, key.sigma); // error term
    for (int i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorus32Distrib(rng);
        tlweSample.b += key.s[i] * tlweSample.a[i];
    }
}

// mu = b - as
double symDecTlweSample(Tlwe& in, const TlweKey& key) {
    auto n = key.n;
    Torus aXs= 0;
    for (int i = 0; i < n; i++) {
        aXs += in.a[i] * key.s[i];
    }
    return torus32ToDouble(in.b - aXs);
}

void rescaleTlweFromTorus32(ScaledTlwe& output, const Tlwe& input) {
    const int n = input.n;
    const int newMod = output.mod;
    output.b = modSwitchFromTorus32(input.b, newMod);
    for (int i = 0; i < n; i++) {
        output.a[i] = modSwitchFromTorus32(input.a[i], newMod);
    }
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
