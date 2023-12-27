//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "tlwe.h"
#include "numeric_functions.h"

using namespace std;

void tlweInitKey(TlweKey& key, const int n, const double sigma) {
    key.n = n;
    key.sigma = sigma;
    key.s = new Integer;
}

void tlweNewBinaryKey(TlweKey& key, const int n, const double sigma) {
    tlweInitKey(key, n, sigma);
    lweKeyGen(key, n);
}

void lweKeyGen(TlweKey& key, const int n) {
    uniform_int_distribution<int64_t> distribution(0, 1);
    std::cout <<"TlweKey: ";
    for (int i = 0; i < n; i++) {
        key.s[i] = distribution(generator);
        std::cout <<key.s[i]<<" ";
    }
    std::cout <<endl;
}

void deleteLweKey(TlweKey& key) {
    delete key.s;
    key.s = nullptr;
}