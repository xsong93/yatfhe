//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/ntt.h"
#include "yatfhe/gadget_decomposition.h"

void trgswFunctionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk, inputModN2, param); //todo:debug
    extractTlweFromTrlwe(tmp, accum, 0); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    tlweKeySwitch(out, ksk, tmp, param);
}

/**
 * Multiply the accumulator by X^sum(bara_i * s_i)
 * */
void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe temp {param.k, param.N};
        controlMux(temp, accum, input.a[i], bsk.bskDft[i], param); // todo:debug
        swap(accum, temp); // assign the previous result to accumulator
    }
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMux(Trlwe& res, const Trlwe& input, const int aBarI, const TrgswDft& bskI, const YatfheParameters& param) {
    trlweRotateMinusOne(res, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    accMulToBsk(res, bskI, param); // res *= bskI // todo: debug
    trlweAccumulate(res, input); // res += input
}

// accum -(GD)> decomp -(ntt)> decompDft -(mul)> accDft -(intt)> accum
void accMulToBsk(Trlwe& accum, const TrgswDft& bskI, const YatfheParameters& param) {
    const auto k = param.k;
    const auto l = param.l;
    const auto N = param.N;
    TrlweDft accDft {k, N};
    DecomposedTrlwe decomp {param};

    gadgetDecomposeTrlwe(decomp, accum, param); // gadget decomposition, G^-1 * TGLWE, T_(N,q)^(k+1) -> Z_N^(k+1)*l

    // ntt
    for (auto lvl = 0; lvl < l; lvl++) {
        applyNttForAB(decomp.rlweDfts[lvl], decomp.rlwes[lvl]);
    }

    // accum += bsk (*) accum, point-wisely
    // https://www.zama.ai/post/tfhe-deep-dive-part-3
    // <Decomp(B), C_k> + Σ_0^(k-1)<Decomp(A_i), C_i>
    // BSK_lrc (*) D_lr = R_c
    for (auto row = 0; row < k + 1; row++) {
        auto& currRes = (row < k) ? accDft.a[row] : accDft.b;
        for (auto lvl = 0; lvl < l; lvl++) {
            for (auto col = 0; col < k; col++) {
                modularAccumulate(currRes.coeffs, decomp.rlweDfts[lvl].a[col].coeffs, bskI.trlweDftSamples[lvl][row].a[col].coeffs);
            }
            modularAccumulate(currRes.coeffs, decomp.rlweDfts[lvl].b.coeffs, bskI.trlweDftSamples[lvl][row].b.coeffs);
        }
    }

    applyInttForAB(accum, accDft); // intt
}

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey) {
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

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey) {
    for (auto i = 0; i < bsk.n; i++) {
        trgswEncrypt(bsk.bsk[i], bsk.bskDft[i], param, trgswKey, tlweKey.s[i]);
    }
}














