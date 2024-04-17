//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"
#include "yatfhe/ntt.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"

// trgsw(0): [trlwe(0)]  (k+1)l rows
void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
//            Trlwe& trlweSample = trgsw.trlweSamples[lvl][row];
//            TrlweDft& trlweDftSample = trgswDft.trlweDftSamples[lvl][row];
//            initCoeffsWithGaussianNoiseSingleSample(trlweSample.b.coeffs, 0, param.lweStdDev); // init b = 0 + e
//            applyNtt(trlweDftSample.b, trlweSample.b);
//            for (auto col = 0; col < param.k; col++) {
//                initCoeffsViaUniformDistribution(trlweSample.a[col].coeffs); // init a
//                applyNtt(trlweDftSample.a[col], trlweSample.a[col]);
//                modularAccumulate(trlweDftSample.b.coeffs, trlweDftSample.a[col].coeffs, trgswKey.trlweKey.sDft[col].coeffs);
//            }
//            applyIntt(trlweSample.b, trlweDftSample.b);
            symEncTrlweSingleSample(trgsw.trlweSamples[lvl][row], trgswDft.trlweDftSamples[lvl][row], trgswKey.trlweKey, 0, param.lweStdDev);
        }
    }
}

// output += mu * G^T
void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
    // add the diagonal matrix (mu * G^T)_ijk to the output
    //       ( 1/B^l                         )
    //      .                              . .
    //    .                              .   .
    //  ( 1/B^2                        )     .
    // ( 1/B                         )       .
    // (     1/B                     )       .
    // (          .                  )       .
    // (              .              )     .
    // (                  .          )   .
    // (                      .      ) .
    // (                         1/B )
    // ( a_0  a_1          a_k-1  b  )

    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits -  (lvl + 1) * param.radixBits);
        // todo: decompose on second level
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
//            LagrangePolynomial tmp(param.N);
//            generateLagrangePolynomialWithValueAt(tmp, mu * g[layer], 0);
//            lagrangePolynomialAccumulate(output.trlweDftSamples[layer][row].b, tmp);
        }
    }
}

//// trgsw(0): [trlwe(0)]  (k+1)l rows
//void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey) {
//    for (auto lvl = 0; lvl < param.lDft; lvl++) {
//        for (auto row = 0; row < param.k + 1; row++) {
//            symEncTrlweSingleSample(trgsw.trlweSamples[lvl][row], trgswDft.trlweDftSamples[lvl][row], trgswKey.trlweKey, 0, param.lweStdDev);
//        }
//    }
//}
//
//// output += mu * G^T
//void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, int64_t mu, const YatfheParameters& param) {
//    // add the diagonal matrix (mu * G^T)_ijk to the output
//    //       ( 1/B^l                         )
//    //      .                              . .
//    //    .                              .   .
//    //  ( 1/B^2                        )     .
//    // ( 1/B                         )       .
//    // (     1/B                     )       .
//    // (          .                  )       .
//    // (              .              )     .
//    // (                  .          )   .
//    // (                      .      ) .
//    // (                         1/B )
//    // ( a_0  a_1          a_k-1  b  )
//
//    for (auto lvl = 0; lvl < param.lDft; lvl++) {
//        auto decomposedMu = mu << (param.dftBits -  (lvl + 1) * param.radixBits);
//        cout << decomposedMu << endl;
//        LagrangePolynomial tmp(param.N, decomposedMu);
//        for (auto row = 0; row < param.k + 1; row++) {
//
//            // add to a_lii
//            if (row < param.k) {
//                lagrangePolynomialAdd(trgswDft.trlweDftSamples[lvl][row].a[row], trgswDft.trlweDftSamples[lvl][row].a[row], tmp);
//                applyIntt(trgsw.trlweSamples[lvl][row].a[row], trgswDft.trlweDftSamples[lvl][row].a[row]);
//                continue;
//            }
//
//            // add to b_lk
//            lagrangePolynomialAdd(trgswDft.trlweDftSamples[lvl][row].b, trgswDft.trlweDftSamples[lvl][row].b, tmp);
//            applyIntt(trgsw.trlweSamples[lvl][row].b, trgswDft.trlweDftSamples[lvl][row].b);
//        }
//    }
//}

/**
 * To encrypt a message as a trgsw ct, there are two steps: 1) generate a trgsw ct with each level and row a rlwe
 * encryption of zero. 2) add mu * G^T to the above trgsw ct
 * @param trgsw Trgsw encryption of message.
 * @param trgswDft Trgsw encryption of message under ntt domain.
 * @param param YatfheParameters
 * @param trgswKey TrgswKey
 * @param mu Message.
 */
void trgswEncrypt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, TrgswKey& trgswKey, const Integer mu) {
    trgswEncZeroNtt(trgsw, trgswDft, param, trgswKey);
    trgswAddIntegerNtt(trgswDft, trgsw, mu, param);
}

// To decrypt, it is sufficient to decrypt the last GLev ciphertext, which is a GLev encryption of m/B^l.
Integer trgswDecrypt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRounding(tmp, trgswDft.trlweDftSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.lweStdDev) >> (param.torusBits - param.radixBits);
}

void trgswExternalProduct(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto l = trgswInput.l;
    DecomposedTrlwe decomposedTrlwe {param};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    // todo: fix overflow issue
    for (auto row = 0; row < k + 1; row++) {
        for (auto lvl = 0; lvl < l; lvl++) {
            for (auto col = 0; col < k; col++) {
                polynomialMulAccNaive(output.a[col], decomposedTrlwe.rlwes[lvl].a[col], trgswInput.trlweSamples[lvl][row].a[col]);
            }
            polynomialMulAccNaive(output.b, decomposedTrlwe.rlwes[lvl].b, trgswInput.trlweSamples[lvl][row].b);
        }
    }
}

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto l = trgswInput.l;
    const auto N = trlweInput.b.N;
    TrlweDft trlweDft {k, N};
    DecomposedTrlwe decomposedTrlwe {param};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    // ntt
    for (auto lvl = 0; lvl < l; lvl++) {
        applyNttForAB(decomposedTrlwe.rlweDfts[lvl], decomposedTrlwe.rlwes[lvl]);
        printArray(decomposedTrlwe.rlwes[lvl].b.coeffs, "b0");
//        printArray(decomposedTrlwe.rlweDfts[lvl].b.coeffs, "bNtt");
        applyInttForAB(decomposedTrlwe.rlwes[lvl], decomposedTrlwe.rlweDfts[lvl]);
        printArray(decomposedTrlwe.rlwes[lvl].b.coeffs, "b1");
    }
//    // accum += bsk (*) accum, point-wisely
//    // https://www.zama.ai/post/tfhe-deep-dive-part-3
//    // <Decomp(B), C_k> + Σ_0^(k-1)<Decomp(A_i), C_i>
//    // BSK_lrc (*) D_lr = R_c
//    for (auto lvl = 0; lvl < l; lvl++) {
//        for (auto row = 0; row < k + 1; row++) {
//            auto& currRes = (row < k) ? decomposedTrlwe.rlweDfts[lvl].a[row] : decomposedTrlwe.rlweDfts[lvl].b;
//            for (auto col = 0; col < k; col++) {
//                modularAccumulate(currRes.coeffs, decomposedTrlwe.rlweDfts[lvl].a[col].coeffs, trgswInput.trlweDftSamples[lvl][row].a[col].coeffs);
//            }
//            modularAccumulate(currRes.coeffs, decomposedTrlwe.rlweDfts[lvl].b.coeffs, trgswInput.trlweDftSamples[lvl][row].b.coeffs);
//        }
//    }
//    // intt
//    for (auto lvl = 0; lvl < l; lvl++) {
//        applyInttForAB(decomposedTrlwe.rlwes[lvl], decomposedTrlwe.rlweDfts[lvl]);
//    }
    recomposeTrlwe(output, decomposedTrlwe, param);



//    applyNttForAB(trlweDft, trlweInput);
//    gadgetDecomposeTrlweNtt(decomposedTrlwe, trlweDft, param);
//
//    //    // accum += bsk (*) accum, point-wisely
////    // https://www.zama.ai/post/tfhe-deep-dive-part-3
////    // <Decomp(B), C_k> + Σ_0^(k-1)<Decomp(A_i), C_i>
////    // BSK_lrc (*) D_lr = R_c
////    for (auto lvl = 0; lvl < l; lvl++) {
////        for (auto row = 0; row < k + 1; row++) {
////            auto& currRes = (row < k) ? decomposedTrlwe.rlweDfts[lvl].a[row] : decomposedTrlwe.rlweDfts[lvl].b;
////            for (auto col = 0; col < k; col++) {
////                modularAccumulate(currRes.coeffs, decomposedTrlwe.rlweDfts[lvl].a[col].coeffs, trgswInput.trlweDftSamples[lvl][row].a[col].coeffs);
////            }
////            modularAccumulate(currRes.coeffs, decomposedTrlwe.rlweDfts[lvl].b.coeffs, trgswInput.trlweDftSamples[lvl][row].b.coeffs);
////        }
////    }
//
//    recomposeTrlweNtt(trlweDft, decomposedTrlwe, param);
//    applyInttForAB(output, trlweDft);
}


