//
// Created by Xintong Song on 2024/3/8.
//
#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/keyswitching.h"

/**
 * The key switching key is a LWE encryption of decomposed current secret key under target secret key.
 * ksk = ksk[i,j] <- TLWE_s'(s_i * B^-j), where i in [0, j*N), j in [0, ks_level)
 * @param ksk TlweKeySwitchingKey
 * @param currKey TrlweKey: Current key to be switched.
 * @param targetKey TlweKey: Target key remaining afterwards.
 * @param param YatfheParameters
 */
void genTlweKeySwitchingKey(TlweKeySwitchingKey& ksk, const TrlweKey& currKey, const TlweKey& targetKey, const YatfheParameters& param) {
    TlweKey inKey(param.k * param.N);
    convertTrlweKeyToTlweKey(inKey, currKey);
    for (auto i = 0; i < inKey.n; i++) {
        auto sOverB = decomposeOverB(inKey.s[i], param); // s_i * B^-j
        for (auto j = 0; j < param.ksLevel; j++) {
            symEncTlweSample(ksk.decomposedKsk[i][j], sOverB[j], targetKey); // encrypt decomposed s under target secret key, no need to map it to Torus, since it's either 0 or 1
        }
    }
}

/**
 * Basically, the idea is to homomorphically cancel the current secret key and re-encrypt it under a new secret key.
 * @param output Tlwe: Tlwe ciphertext under new key.
 * @param ksk TlweKeySwitchingKey
 * @param input Tlwe: Tlwe ciphertext to be keyswitched.
 * @param param YatfheParameters
 */
void tlweKeySwitch(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param) {
    output.b = input.b; // init output as (0,..., 0, b)
    for (auto i = 0; i < input.n; i++) {
        vector<Torus> aBar(param.ksLevel);
        Tlwe tmp(output.n);
        signedGadgetDecomposition(aBar, input.a[i], param); // todo: (aBar_1, ..., aBar_l) <- g^-1(ai)
        for (auto j = 1; j <= param.ksLevel; j++) {
            // todo: dot(aj, kskij)
            for (auto k = 0; k < input.n; k++) {
                tmp.a[k] += aBar[j] * ksk.decomposedKsk[i][j].a[k];
            }
            tmp.b += aBar[j] * ksk.decomposedKsk[i][j].b;
        }
        lweSubTo(output, tmp);
    }
}