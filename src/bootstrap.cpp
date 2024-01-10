//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "bootstrap.h"
#include "numeric_functions.h"
#include "fft.h"

void newBootstrappingKey(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey, const int unfolding) {
    if (unfolding == 1) {
        newBootstrappingKeyWoUnfolding(bsk, trgswKey, tlweKey);
    }
////    const int l = trgswKey->l, Bg_bit = trgswKey->Bg_bit, k = trgswKey->trlwe_key->k, N = trgswKey->trlwe_key->bskDft[0]->N;
//    BootstrappingKey* res{new BootstrappingKey};
////    res->n = tlweKey->n;
////    res->k = k;
////    res->l = l;
////    res->N = N;
////    res->Bg_bit = Bg_bit;
//    res->unfolding = unfolding;
//    const int keyExp = 1 << unfolding;
//    const int finalExp = keyExp / unfolding; // expansion constants
//    res->su = trgsw_alloc_new_sample_array(tlweKey->n*finalExp, l, Bg_bit, k, N);
//
//    for (int i = 0; i < tlweKey->n; i+=unfolding){
//        for (int j = 0; j < keyExp; j++){
//            Binary key = 1;
//            for (int u = 0, j_ = j; u < unfolding; u++, j_>>=1){
//                if(j_&1) key *= tlweKey->s[i + u];
//                else     key *= 1 - tlweKey->s[i + u];
//            }
//            trgsw_monomial_sample(res->su[i*finalExp + j], key, 0, trgswKey);
//        }
//    }
}

void newBootstrappingKeyWoUnfolding(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey) {
    const int n = tlweKey.n;
    bsk.bskDft.resize(n);
    bsk.unfolding = 1;
    bsk.bsk.resize(n);
    for (int i = 0; i < n; i++) {
        Trgsw* trgsw = &bsk.bsk[i];
        initTrgswDftSample(bsk.bskDft[i], trgswKey);
        initTrgswSample(*trgsw, trgswKey);
        trgswEncZero(*trgsw, tlweKey.sigma, trgswKey);
        const Integer message = tlweKey.s[i];
        //tGswAddMuIntH(result, message, key->params);
    }
}

// trgsw(0)
void trgswEncZero(Trgsw& trgsw, const double sigma, const TrgswKey& trgswKey) {
    const int N = trgswKey.trlweKey.s[0].N;
    const int k = trgswKey.trlweKey.k;
    const int l = trgswKey.l;
    const int kpl = (k + 1) * l;
    for (int p = 0; p < kpl; p++) {
        Trlwe* trlwe = &trgsw.trlweSamples[p];
        for (int j = 0; j < N; j++) {
            trlwe->b.coeffs[j] = addGaussianNoise(0, sigma);
        }
//        fft(trlwe.b);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < N; j++) {
                trlwe->a[i].coeffs[j] = uniformTorus32Distrib(rng);
            }
//            torusPolynomialAddMulR(trlwe.b, trgswKey.trlweKey.s[i], trlwe.a[i]);
        }
    }
}











