//
// Created by Xintong Song on 2024/3/8.
//
#include "yatfhe/gadget_decomposition.h"
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
    TlweKey inKey {param.k * param.N};
    convertTrlweKeyToTlweKey(inKey, currKey);
    for (auto i = 0; i < inKey.n; i++) {
        std::vector<Integer> sOverB(param.ksLevel);
        decomposeOverB(sOverB, inKey.s[i], param); // s_i * B^-j
        for (auto j = 0; j < param.ksLevel; j++) {
            symEncTlwe(ksk.decomposedKsk[i][j], sOverB[j], targetKey);
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
void switchKeyForTlwe(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param) {
    resetTlweToZero(output);
    output.b = input.b; // init output as (0,..., 0, b)
    Torus mul = 0;
    for (auto i = 0; i < input.n; i++) {
        DecomposedData aBar {param.ksLevel};
        Tlwe tmp {output.n};
//        signedGadgetDecomposition(aBar, input.a[i], param); // (aBar_1, ..., aBar_l) <- g^-1(ai)
        gadgetDecompose(aBar, input.a[i], param);
        for (auto j = 0; j < param.ksLevel; j++) {
            for (auto k = 0; k < output.n; k++) {
                mul = multTorus(LWE_Q, aBar.value[j] * aBar.sign, ksk.decomposedKsk[i][j].a[k]);
                tmp.a[k] = addTorus(LWE_Q, tmp.a[k], mul);
            }
            mul = multTorus(LWE_Q, aBar.value[j] * aBar.sign, ksk.decomposedKsk[i][j].b);
            tmp.b = addTorus(LWE_Q, tmp.b, mul);
        }
        subTlweInPlace(output, tmp);
    }
}