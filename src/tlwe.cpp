//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yatfhe_parameters.h"
#include "tlwe.h"
#include "numeric_functions.h"

using namespace std;

void initTlweKey(TlweKey& key, const int n, const double sigma) {
    key.n = n;
    key.sigma = sigma;
//    key.bskDft = new Integer[n];
    key.s.resize(n);
}

void newBinaryTlweKey(TlweKey& key, const YatfheParameters& parameters) {
    initTlweKey(key, parameters.n, parameters.lweStdDev);
    lweKeyGen(key, parameters.n);
}

void lweKeyGen(TlweKey& key, const int n) {
    uniform_int_distribution<int> distribution(0, 1);
    std::cout <<"TlweKey: ";
    for (int i = 0; i < n; i++) {
        key.s[i] = distribution(rng);
        std::cout << i << ":" <<key.s[i]<<" ";
    }
    std::cout <<endl;
}

void deleteLweKey(TlweKey& key) {
//    delete key.bskDft;
//    key.bskDft = nullptr;
}