//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yautil/numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

void lweKeyGen(TlweKey& key, const int n) {
    uniform_int_distribution<int> distribution(0, 1);
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

// Basically, the idea is to homomorphically cancel the secret key and re-encrypt it under a new secret key.
void lweKeySwitch(Tlwe& output, const Tlwe& keySwitchingKey, Tlwe& input, const YatfheParameters& param) {
    output.b = input.b; // init output as (0,..., 0, b)
    auto precOffset = 1 << (param.torusBits - (1 + param.radixBits * param.t)); //precision
    auto g = genGadgetVector(param.radixBits, param.t, param.torusBits);
    for (auto i = 0; i < param.n; i++) {
        uint32_t barai = input.a[i] + precOffset;

        // signed decomp
        for (auto j = 1; j <= param.t; j++) {
            auto aij = (barai >> (param.torusBits - j * param.radixBits)) & param.digitMask;
            if (aij != 0) {
                lweSubTo(output, keySwitchingKey);
            }
        }
    }
}

// output -= input
void lweSubTo(Tlwe& output, const Tlwe& input) {
    auto n = output.n;
    output.b -= input.b;
    for (auto i = 0; i < n; i++) {
        output.a[i] -= input.a[i];
    }
}

void tlweCopy(Tlwe& output, const Tlwe& input) {
    const auto n = input.n;
    for (auto i = 0; i < n; i++) {
        output.a[i] = input.a[i];
    }
    output.b = input.b;
}
