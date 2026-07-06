//
// Created by Xintong Song on 2024/3/8.
//
#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
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

void genPrivateKeySwitchingKey(PrivateKeySwitchingKey& psk, const TrlweKey& trlweKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    const auto n = param.n;
    // f(m) is treated as a monomial in R; embed scalar as X^0.
    const int monoIndex = 0;
    encTrlevSingleSampleMonomial(psk.one, trlweKey, 1, monoIndex, param);
    for (auto i = 0; i < n; i++) {
        encTrlevSingleSampleMonomial(psk.trlevs[i], trlweKey, tlweKey.s[i], monoIndex, param);
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
    const int nOut = output.n;

    // Accumulate all KSK contributions in int64 with no per-element modular reduction.
    // Max magnitude: input.n * ksLevel * (radixBase/2) * (LWE_Q/2) <= 2048*4*128*2^19 = 2^40, fits in int64.
    vector<int64_t> delta(nOut, 0);
    int64_t delta_b = 0;

    DecomposedData aBar{param.ksLevel};  // hoisted: avoids 2*input.n heap allocations in the loop
    for (int i = 0; i < input.n; i++) {
        gadgetDecompose(aBar, input.a[i], param);
        for (int j = 0; j < param.ksLevel; j++) {
            const int64_t coeff = static_cast<int64_t>(aBar.value[j]) * aBar.sign;
            if (coeff == 0) {
                continue;
            }
            const auto& kskij = ksk.decomposedKsk[i][j];  // hoist double dereference
            for (int k = 0; k < nOut; k++) {
                delta[k] += coeff * kskij.a[k];
            }
            delta_b += coeff * kskij.b;
        }
    }

    // Apply with one reduction per element.
    for (int k = 0; k < nOut; k++) {
        output.a[k] = static_cast<Torus>(longModP(-delta[k], LWE_Q));
    }
    output.b = static_cast<Torus>(longModP(static_cast<int64_t>(input.b) - delta_b, LWE_Q));
}

// LWE-to-RLWE Private Key Switching (Wang et al., 2024, Sec. 2.3.1).
// PrivateKS_f(c=(a,b)) = sum_i a_i ⊙ RLWE'(f(s_i)) + b ⊙ RLWE'(f(1)).
// Here f is identity and RLWE'(f(s_i)) is reused from bootstrap key cPrime terms.
void tlweToTrlwePrivateKeySwitching(Trlwe& out, const Tlwe& in, const PrivateKeySwitchingKey& psk,
                                    const YatfheParameters& param) {
    // Ensure output is clean; multTrlevWithConst adds into output.
    clearTrlwe(out);

    // b ⊙ RLWE'(1): represented as a noiseless RLWE with only constant term.
    // `in` is in LWE torus domain, `out` is in RLWE torus domain -> rescale.
    const auto bScaled = modSwitchFromTorusGeneral(in.b, TORUS_Q, LWE_Q);
    multTrlevWithConst(out, psk.one, bScaled, param);

    // -sum_i a_i ⊙ RLWE'(s_i): TLWE phase is b - <a, s>, so the ai terms must be subtracted.

    for (int i = 0; i < in.n; ++i) {
        if (in.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        auto& currKey = psk.trlevs[i];
        const Torus a_i_scaled = modSwitchFromTorusGeneral(in.a[i], TORUS_Q, LWE_Q);
        multTrlevWithConst(tmp, currKey, a_i_scaled, param);
        subTrlwe(out, out, tmp);
    }
}