//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "bootstrap.h"
#include "numeric_functions.h"
#include "ntt.h"


void trgswFunctionalBootstrapping(TrgswDft& out, const Tlwe& in, const BootstrappingKey& bsk, const Torus msg, const YatfheParameters& param) {
    const int n = param.n;
    const int N = param.N;
    const int l = param.l;
    const int bgBit = param.bgBit;
    const int k = param.k;
    const int N2 = N * 2;
    const int logN2 = (int) log2(N2);
    const int torusBase = param.torusBase;
    const Torus precOffset = doubleToTorus32(1.0 / (4 * torusBase));
    NegaCyclicTlwe negaCyclicInput(N2);
    modSwitchFromTorus32ToN2(negaCyclicInput, in);
    Trgsw tv(param);
    Trlwe acc(k + 1, N);
//    genNoiselessTrgswSample(tv, msg, param);
    genNoiselessTrlweSample(acc, msg, negaCyclicInput, param);
    Trgsw tmp(param);
//    trgsw_mul_bly_xai(tmp, tv, N2 - torus2int(in->b + precOffset, logN2));
    blindRotate(acc, bsk, negaCyclicInput, param);
//    trgsw_to_DFT(out, tmp);
}

void blindRotate(Trlwe& acc, const BootstrappingKey& bsk, const NegaCyclicTlwe& sample, const YatfheParameters& param) {
    const int n = param.n;
    const int k = param.k;
    const int N = param.N;
    const auto bara = sample.a;
    Trlwe temp(k + 1, N);
    for (int i = 0; i < n; i++) {
        if (bara[i] == 0) {
            continue;
        }
        muxRotate(temp, acc, bsk.bskDft[i], bara[i], param);
    }

    // todo
}

// ACC = BSKi * [(X^barai - 1) * ACC] + ACC
void muxRotate(Trlwe& res, Trlwe& acc, const TrgswDft& bski, const int barai, const YatfheParameters& param) {
    const auto k = param.k;

    // res = (X^barai - 1) * ACC
    for (int i = 0; i <= k; i++) {
        torusPolynomialMulByXaiMinusOne(res.a[i], barai, acc.a[i]);
    }

    // acc *= BKi
    trgswMulToTrlwe(acc, bski, param);

    // todo
}

// accum -(GD)> deca -(fft)> decaFFT -(mul)> tmpa -(ifft)> accum
void trgswMulToTrlwe(Trlwe& acc, const TrgswDft& bski, const YatfheParameters& param) {
    const int k = param.k;
    const int l = param.l;
    const int N = param.N;
    const int kpl = (k + 1) * l;
    vector<vector<IntPolynomial>> decomp(k + 1, vector<IntPolynomial>(l, IntPolynomial(N)));
    vector<vector<LagrangePolynomial>> decompDft(k + 1, vector<LagrangePolynomial>(l, LagrangePolynomial(N)));

    // gadget decomposition, G^-1 * TGLWE, T_(N,q)^(k+1) -> Z_N^(k+1)*l
    gadgetDecomposition(decomp , acc.a, param);

    // ntt
    for (int i = 0; i <= k; i++) {
        for (int p = 0; p < l; p++) {
            applyNtt(decompDft[i][p], decomp[i][p]);
        }
    }

    // acc = gsw (*) acc, point-wisely
    for (int i = 0; i <= k; i++) {
        for (int p = 0; p < l; p++) {
            // todo
        }
    }
}

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey, const TlweKey& tlweKey) {
    if (bsk.unfolding == 1) {
        bootstrappingKeyGenWoUnfolding(bsk, param, trgswKey, tlweKey);
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

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey, const TlweKey& tlweKey) {
    const int n = bsk.n;
    for (int i = 0; i < n; i++) {
        Trgsw& trgsw = bsk.bsk[i];
        TrgswDft& trgswDft = bsk.bskDft[i];
//        initTrgswDftSample(trgswDft, param);
//        initTrgswSample(trgsw, param);
        trgswEncZero(trgsw, trgswDft, param, trgswKey);
        // const Integer message = tlweKey.s[i];
        //tGswAddMuIntH(result //trgsw, message, key->params);
    }
}

// trgsw(0)
void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const int N = param.N;
    const int k = param.k;
    const int l = param.l;
    const int kpl = (k + 1) * l;
    const double sigma = param.lweStdDev;
    for (int p = 0; p < kpl; p++) {
        Trlwe& trlweSample = trgsw.trlweSamples[p];
        TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[p];
        initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, N, sigma);
        applyNtt(trlweDftSample.b, trlweSample.b);
        calModularInnerProduct(trlweDftSample.b, trlweSample.a, trgswKey.trlweKey.s, N, k);
    }
}












