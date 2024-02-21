//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "trlwe.h"
#include "bootstrapping.h"
#include "numeric_functions.h"
#include "ntt.h"
#include "yautil/tool.h"
#include "yautil/time_counter.h"

void trgswFunctionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TorusPolynomial& v, const YatfheParameters& param) {
    const int n = param.n;
    const int N = param.N;
    const int k = param.k;
    const int N2 = N * 2;
    ScaledTlwe inputModN2(N2, n);
    Trlwe accum(k + 1, N);
    rescaleTlweFromTorus32(inputModN2, input); //rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk, inputModN2, param);
    extractTlweFromTrlwe(out, accum, 0); // out = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    // todo: keyswitching
//    printTrlweAB(accum, "accum");
}

/**
 * Multiply the accumulator by X^sum(bara_i * s_i)
 * */
void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    const int n = param.n;
    const int k = param.k;
    const int N = param.N;
    const auto bara = input.a;
    for (int i = 0; i < n; i++) {
        if (bara[i] == 0) {
            continue;
        }
        Trlwe temp(k + 1, N);
        trlweRotateMinusOne(temp, accum, bara[i]); // temp = c1 - c0 = (X^barai - 1) * input
        controlMux(temp, accum, bsk.bskDft[i], param);
        swap(temp.a, accum.a);
    }
}

// res = bsk * (c1 - c0) + c0 = bski * [(X^barai - 1) * input] + input
void controlMux(Trlwe& res, const Trlwe& input, const TrgswDft& bski, const YatfheParameters& param) {
    const int k = param.k;

    // res *= bski
    accMulToBsk(res, bski, param); // todo: debug

    // res += input
    for (int i = 0; i < k + 1; i++) {
        polynomialAccumulate(res.a[i], input.a[i]);
    }
}

// accum -(GD)> decomp -(ntt)> decompDft -(mul)> accDft -(intt)> accum
void accMulToBsk(Trlwe& accum, const TrgswDft& bski, const YatfheParameters& param) {
    const int k = param.k;
    const int l = param.l;
    const int N = param.N;
    vector<LagrangePolynomial> accDft(k + 1, LagrangePolynomial(N));
    vector<vector<IntPolynomial>> decomp(k + 1, vector<IntPolynomial>(l, IntPolynomial(N)));
    vector<vector<LagrangePolynomial>> decompDft(k + 1, vector<LagrangePolynomial>(l, LagrangePolynomial(N)));

    // gadget decomposition, G^-1 * TGLWE, T_(N,q)^(k+1) -> Z_N^(k+1)*l
    gadgetDecomposition(decomp , accum.a, param);
//    printPolyMat(decomp, "decomp");

    // ntt
    for (int i = 0; i < k + 1; i++) {
        for (int p = 0; p < l; p++) {
            applyNtt(decompDft[i][p], decomp[i][p]);
        }
    }

    // accum += bsk (*) accum, point-wisely
    for (int i = 0; i < k + 1; i++) {
        for (int j = 0; j < l; j++) {
            for (int m = 0; m < k + 1; m++) {
                modularAccumulate(accDft[m].coeffs, decompDft[i][j].coeffs, bski.trlweDftSamples[i][j].a[m].coeffs);
            }
        }
    }

    // intt
    for (int i = 0; i < k + 1; i++) {
        applyIntt(accum.a[i], accDft[i]);
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
    for (int i = 0; i < k + 1; i++) {
        for (int j = 0; j < l; j++) {
            Trlwe& trlweSample = trgsw.trlweSamples[i][j];
            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[i][j];
            initCoeffsWithGaussianNoise(trlweSample.b.coeffs, 0, N, sigma); // init b
            applyNtt(trlweDftSample.b, trlweSample.b);
            for (int m = 0; m < k; m++) {
                initCoeffsViaUniformDistribution(trlweSample.a[m].coeffs, N); // init a
                applyNtt(trlweDftSample.a[m], trlweSample.a[m]);
                calModularInnerProductNtt(trlweDftSample.b, trlweDftSample.a[m], trgswKey.trlweKey.s[m], N);
            }
        }
    }
}












