//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe_parameters.h"
#include "tlwe.h"
#include "numeric_functions.h"
#include "yautil/tool.h"

using namespace std;

//void initTlweKey(TlweKey& key, const int n, const double sigma) {
//    key.n = n;
//    key.sigma = sigma;
////    key.bskDft = new Integer[n];
//    key.s.resize(n);
//}

//void newBinaryTlweKey(TlweKey& key, const YatfheParameters& param) {
//    initTlweKey(key, param.n, param.lweStdDev);
//    lweKeyGen(key, param.n);
//}

void lweKeyGen(TlweKey& key, const int n) {
    uniform_int_distribution<int> distribution(0, 1);
    for (int i = 0; i < n; i++) {
        key.s[i] = distribution(rng);
    }
//    printArray(key.s, "TlweKey");
}

//void initTlweSample(Tlwe& tlwe, int n) {
//    tlwe.a.resize(n);
//    tlwe.n = n;
//}

// b = aj * sj + u + e
void symEncTlweSample(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    tlweSample.b = addGaussianNoise(message, key.sigma); // error term
    for (int i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorus32Distrib(rng);
        tlweSample.b += key.s[i] * tlweSample.a[i];
    }
}

// mu = b - as
double symDecTlweSample(Tlwe& in, TlweKey& key) {
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

//void deleteLweKey(TlweKey& key) {
////    delete key.bskDft;
////    key.bskDft = nullptr;
//}