//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "bootstrap.h"
#include "numeric_functions.h"

BootstrappingKey* newBootstrappingKey(const TrgswKey* trgswKey, const TlweKey* tlweKey, const int unfolding) {
    if (unfolding == 1) {
        return newBootstrappingKeyWoUnfolding(trgswKey, tlweKey);
    }
//    const int l = trgswKey->l, Bg_bit = trgswKey->Bg_bit, k = trgswKey->trlwe_key->k, N = trgswKey->trlwe_key->s[0]->N;
    BootstrappingKey* res{new BootstrappingKey};
//    res->n = tlweKey->n;
//    res->k = k;
//    res->l = l;
//    res->N = N;
//    res->Bg_bit = Bg_bit;
    res->unfolding = unfolding;
    const int keyExp = 1 << unfolding;
    const int finalExp = keyExp / unfolding; // expansion constants
    res->su = trgsw_alloc_new_sample_array(tlweKey->n*finalExp, l, Bg_bit, k, N);

    for (int i = 0; i < tlweKey->n; i+=unfolding){
        for (int j = 0; j < keyExp; j++){
            Binary key = 1;
            for (int u = 0, j_ = j; u < unfolding; u++, j_>>=1){
                if(j_&1) key *= tlweKey->s[i + u];
                else     key *= 1 - tlweKey->s[i + u];
            }
            trgsw_monomial_sample(res->su[i*finalExp + j], key, 0, trgswKey);
        }
    }
    return res;
}

BootstrappingKey* newBootstrappingKeyWoUnfolding(const TrgswKey* trgswKey, const TlweKey* tlweKey) {
    BootstrappingKey* bsk{new BootstrappingKey};
    bsk->s = new TrgswDft;
    initTrgswDftSample(bsk->s, tlweKey->n, trgswKey->l, trgswKey->bgBit, trgswKey->trlweKey->k, trgswKey->trlweKey->s[0].N);
    bsk->unfolding = 1;
    Trgsw* trgsw{new Trgsw};
    initTrgswSample(trgsw, trgswKey->l, trgswKey->bgBit, trgswKey->trlweKey->k, trgswKey->trlweKey->s[0].N);
    for (int i = 0; i < tlweKey->n; i++) {
//        trgsw_monomial_sample(trgsw, tlweKey->s[i], 0, trgswKey);
//        trgsw_to_DFT(res->s[i], trgsw);
        trgswEncZero(trgsw, tlweKey->sigma, trgswKey);
    }
    free_trgsw(trgsw);
    return bsk;
}

void trgswEncZero(Trgsw *trgsw, const double sigma, const TrgswKey* key) {
    const int N = key->trlweKey->s[0].N;
    const int k = key->trlweKey->k;
    const int l = key->l;
    for (int p = 0; p < (k + 1) * l; p++) {
        Trlwe* trlwe = &trgsw->samples[p];
        for (int j = 0; j < N; j++) {
            trlwe->b->coeffs[j] = gaussian32(0, sigma);
        }
        for (int i = 0; i < k; i++) {
            torusPolynomialUniform(&trlwe->a[i]);
            torusPolynomialAddMulR(trlwe->b, &key->trlweKey[i], &trlwe->a[i]);
        }
    }
}











