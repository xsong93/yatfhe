//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "bootstrap.h"
#include "numeric_functions.h"
#include "ntt.h"

void newBootstrappingKey(BootstrappingKey& bsk, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey, const TlweKey& tlweKey) {
    if (bsk.unfolding == 1) {
        newBootstrappingKeyWoUnfolding(bsk, yatfheParameters, trgswKey, tlweKey);
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

void newBootstrappingKeyWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey, const TlweKey& tlweKey) {
    const int n = bsk.n;
    for (int i = 0; i < n; i++) {
        Trgsw& trgsw = bsk.bsk[i];
        TrgswDft& trgswDft = bsk.bskDft[i];
//        initTrgswDftSample(trgswDft, yatfheParameters);
//        initTrgswSample(trgsw, yatfheParameters);
        trgswEncZero(trgsw, trgswDft, yatfheParameters, trgswKey);
        // const Integer message = tlweKey.s[i];
        //tGswAddMuIntH(result //trgsw, message, key->params);
    }
}

// trgsw(0)
void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey) {
    const int N = yatfheParameters.N;
    const int k = yatfheParameters.k;
    const int l = yatfheParameters.l;
    const int kpl = (k + 1) * l;
    const double sigma = yatfheParameters.lweStdDev;
    for (int p = 0; p < kpl; p++) {
        Trlwe& trlweSample = trgsw.trlweSamples[p];
        TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[p];
        addGaussianNoiseToCoeffs(trlweSample.b.coeffs, N, sigma);
        applyNtt(trlweSample.b, trlweDftSample.b);
        for (int i = 0; i < k; i++) {
            initCoeffsViaUniformDistribution(trlweSample.a[i].coeffs, N);
            LagrangePolynomial sDft {N};
            LagrangePolynomial aDft {N};
            applyNtt(trgswKey.trlweKey.s[i], sDft);
            applyNtt(trlweSample.a[i], aDft);
            calculateB(aDft.coeffs, sDft.coeffs, trlweDftSample.b.coeffs, N);
        }
    }
}

void trgswFunctionalBootstrapping(TrgswDft& out, const Tlwe& in, const BootstrappingKey& bsk, const YatfheParameters& parameters) {
    const int N = parameters.N;
    const int l = parameters.l;
    const int bgBit = parameters.bgBit;
    const int k = parameters.k;
    const int N2 = N * 2;
    const int logN2 = (int) log2(N2);
    const int torusBase = parameters.torusBase;
    const Torus precOffset = doubleToTorus32(1.0 / (4 * torusBase));
    Trgsw tv = trgsw_new_noiseless_trivial_sample(1, l, bgBit, k, N);
    Trgsw tmp = trgsw_alloc_new_sample(l, bgBit, k, N);
    trgsw_mul_bly_xai(tmp, tv, N2 - torus2int(in->b + precOffset, logN2));
    blind_rotate_trgsw(tmp, in->a, key->s, in->n);
    trgsw_to_DFT(out, tmp);
}

// b = aj*sj
void calculateB(std::vector<uint64_t>& coeffsA, std::vector<uint64_t>& coeffsS, std::vector<uint64_t>& coeffsB, const int N) {
    for (int j = 0; j < N; j++) {
        auto tmp = modMul(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modAdd(coeffsB[j], tmp);
    }
}

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, const int N) {
    for (int j = 0; j < N; j++) {
        coeffs[j] = uniformTorus32Distrib(rng);
    }
}

void addGaussianNoiseToCoeffs(std::vector<Torus>& coeffs, const int N, const double sigma) {
    for (int j = 0; j < N; j++) {
        coeffs[j] = addGaussianNoise(0, sigma);
    }
}











