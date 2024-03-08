//
// Created by Xintong Song on 2024/3/8.
//
#include "yautil/numeric_functions.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/keyswitching.h"

void genTlweKeySwitchingKey(TlweKeySwitchingKey& ksk, const TlweKey& currKey, const TlweKey& newKey, const YatfheParameters& param) {

}

// Basically, the idea is to homomorphically cancel the secret key and re-encrypt it under a new secret key.
void tlweKeySwitch(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param) {
    output.b = input.b; // init output as (0,..., 0, b)
    auto precOffset = 1 << (param.torusBits - (1 + param.radixBits * param.ksLevel)); //precision
    auto g = genGadgetVector(param.radixBits, param.ksLevel, param.torusBits);
    for (auto i = 0; i < input.n; i++) {
        vector<Torus> aBar(param.ksLevel);
//        uint32_t barai = input.a[i] + precOffset;
        signedGadgetDecomposition(aBar, input.a[i]); // todo: g-1(ai)
        Tlwe tmp(output.n);
        // signed decomp
        for (auto j = 1; j <= param.ksLevel; j++) {
//            auto aij = (barai >> (param.torusBits - j * param.radixBits)) & param.digitMask;
            // todo: dot(aij, kskij)
        }
        lweSubTo(output, tmp);
    }
}