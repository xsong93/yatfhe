//
// Created by Xintong Song on 2023/12/25.
//
//#include <omp.h>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/crt.h"

// trgsw(0): [trlwe(0)]  (k+1)l rows
void trgswEncZero(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSample(trgsw.trlweSamples[lvl][row], trgswKey.trlweKey, 0);
        }
    }
}

// output += mu * G^T
void trgswAddInteger(Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
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
                continue;
            }

            // add to b_lk
            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
        }
    }
}

// trgsw(0): [trlwe(0)]  (k+1)l rows
void trgswEncZeroNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSampleNtt(trgsw.trlweSamples[lvl][row], trgswDft.trlweDftSamples[lvl][row], trgswKey.trlweKey, 0);
        }
    }
}

// output += mu * G^T
void trgswAddIntegerNtt(TrgswDft& trgswDft, Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
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
        }
    }
}

void trgswEncrypt(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    trgswEncZero(trgsw, param, trgswKey);
    trgswAddInteger(trgsw, mu, param);
}

/**
 * To encrypt a message as a trgsw ct, there are two steps: 1) generate a trgsw ct with each level and row a rlwe
 * encryption of zero. 2) add mu * G^T to the above trgsw ct
 * @param trgsw Trgsw encryption of message.
 * @param trgswDft Trgsw encryption of message under ntt domain.
 * @param param YatfheParameters
 * @param trgswKey TrgswKey
 * @param mu Message.
 */
void trgswEncryptNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    trgswEncZeroNtt(trgsw, trgswDft, param, trgswKey);
    trgswAddIntegerNtt(trgswDft, trgsw, mu, param);
}

// To decrypt, it is sufficient to decrypt the last GLev ciphertext, which is a GLev encryption of m/B^l.
Integer trgswDecrypt(const Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRounding(tmp, trgsw.trlweSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweStdDev, param.torusBits - param.radixBits);
}

Integer trgswDecryptNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRoundingNtt(tmp, trgswDft.trlweDftSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweStdDev, param.torusBits - param.radixBits);
}

void trgswCRTDecomp(std::vector<TrgswDft24>& out, const Trgsw& in, const YatfheParameters& param) {
    std::vector<Trgsw8> tmp(param.d, Trgsw8{param.dh, param.k, param.N});

    for (size_t d = 0; d < param.d; d++) {
        auto& taoU = param.taoU[d];
        auto& qd = param.qd[d];
        auto& tmpOut = out[d];
        auto& tmpD = tmp[d];
        for (size_t l = 0; l < param.dh; l++) {
            for (size_t k1 = 0; k1 < param.k + 1; k1++) {
                auto& nttOut = tmpOut.trlweDftSamples[l][k1];
                auto& nttIn = tmpD.trlweSamples[l][k1];
                auto& rgswInA = in.trlweSamples[l][k1].a;
                auto& rgswInB = in.trlweSamples[l][k1].b;
                auto& rgswResA = tmpD.trlweSamples[l][k1].a;
                auto& rgswResB = tmpD.trlweSamples[l][k1].b;
                for (size_t k2 = 0; k2 < param.k; k2++) {
                    auto& coeffResA = rgswResA[k2].coeffs;
                    auto& coeffInA = rgswInA[k2].coeffs;
                    for (size_t j = 0; j < param.N; j++) {
                        coeffResA[j] = static_cast<int8_t>(longModP(taoU * coeffInA[j], qd));
                    }
                }
                auto& coeffResB = rgswResB.coeffs;
                auto& coeffInB = rgswInB.coeffs;
                for (size_t j = 0; j < param.N; j++) {
                    coeffResB[j] = static_cast<int8_t>(longModP(taoU * coeffInB[j], qd));
                }
                applyNttForAB24(nttOut, nttIn);
            }
        }
    }
}

void trgswExternalProduct(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level = trgswInput.l;
    DecomposedTrlwe decomposedTrlwe {param};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    // <Decomp(B), C_k^bar> + Σ_(0,k-1)<Decomp(A_i), C_i^bar>
    // https://www.zama.ai/post/tfhe-deep-dive-part-3
    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto& curr = (col < k) ? decomposedTrlwe.rlwes[lvl].a[col] : decomposedTrlwe.rlwes[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                auto& out = (col2 < k) ? output.a[col2] : output.b;
                auto& curr2 = (col2 < k) ? trgswInput.trlweSamples[lvl][col].a[col2]
                                         : trgswInput.trlweSamples[lvl][col].b;
                polynomialMulAccNaive(out, curr, curr2);
            }
        }
    }
}

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswDftInput, Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level = trgswDftInput.l;
    const auto N = trlweInput.b.N;
    TrlweDft trlweDftRes {k, N};
    DecomposedTrlwe decomposedTrlwe {param};
    DecomposedTrlweDft decomposedTrlweDft {param, param.l};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);  // 8 * 2

//#pragma omp parallel for
    for (auto i = 0; i < decomposedTrlwe.l; i++) {
        applyNttForAB(decomposedTrlweDft.rlweDfts[i], decomposedTrlwe.rlwes[i]);
    }

//#pragma omp parallel for collapse(2) private(out)
    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto& curr = (col < k) ? decomposedTrlweDft.rlweDfts[lvl].a[col] : decomposedTrlweDft.rlweDfts[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                auto& out = (col2 < k) ? trlweDftRes.a[col2] : trlweDftRes.b;
                auto& curr2 = (col2 < k) ? trgswDftInput.trlweDftSamples[lvl][col].a[col2]
                                         : trgswDftInput.trlweDftSamples[lvl][col].b;
                calModularInnerProductNtt(out, curr, curr2);
            }
        }
    }
    applyInttForAB(output, trlweDftRes);
}

/*//todo
void trgswExternalProductSplitNtt(Trlwe& output, const TrgswDft& trgswDftInput, Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level = trgswDftInput.l;
    const auto N = trlweInput.b.N;
    TrlweDft trlweDftRes {k, N};
    DecomposedTrlwe decomposedTrlwe {param};
    DecomposedTrlweDft decomposedTrlweDft {param, param.l};

//    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);
    for (auto row = 0; row < k + 1; row++) {
        auto& currIn = (row < k) ? trlweInput.a[row] : trlweInput.b;
        for (auto j = 0; j < N; j++) {
            DecomposedData d {level};
            gadgetDecompose(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < level; lvl++) {
                auto& currOut = (row < k) ? decomposedTrlwe.rlwes[lvl].a[row] : decomposedTrlwe.rlwes[lvl].b;
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }

//#pragma omp parallel for
    for (auto i = 0; i < decomposedTrlwe.l; i++) {
        for (auto row = 0; row < decomposedTrlwe.rlwes[i].a.size(); row++) {
            applyNtt(decomposedTrlweDft.rlweDfts[i].a[row], decomposedTrlwe.rlwes[i].a[row]);
        }
        applyNtt(decomposedTrlweDft.rlweDfts[i].b, decomposedTrlwe.rlwes[i].b);
    }

//#pragma omp parallel for collapse(2) private(out)
    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col = 0; col < k; col++) {
            for (auto col2 = 0; col2 < k + 1; col2++) {
                auto& out = (col2 < k) ? trlweDftRes.a[col2] : trlweDftRes.b;
                auto& curr2 = (col2 < k) ? trgswDftInput.trlweDftSamples[lvl][col].a[col2]
                                         : trgswDftInput.trlweDftSamples[lvl][col].b;
                calModularInnerProductNtt(out, decomposedTrlweDft.rlweDfts[lvl].a[col], curr2);
            }
        }
    }

    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col2 = 0; col2 < k + 1; col2++) {
            auto& out = (col2 < k) ? trlweDftRes.a[col2] : trlweDftRes.b;
            auto& curr2 = (col2 < k) ? trgswDftInput.trlweDftSamples[lvl][k].a[col2]
                                     : trgswDftInput.trlweDftSamples[lvl][k].b;
            calModularInnerProductNtt(out, decomposedTrlweDft.rlweDfts[lvl].b, curr2);
        }
    }
    applyInttForAB(output, trlweDftRes);
}*/

void trgswExternalProductCRT(std::vector<Trlwe8>& output, const std::vector<TrgswDft24>& trgswDftInput, std::vector<Trlwe8>& trlweInput, const YatfheParameters& param) {
    std::vector<Trlwe8> tmpD(param.dh, Trlwe8{param.k, param.N});
    std::vector<std::vector<Trlwe8>> tmpDB(param.d, std::vector<Trlwe8>(param.dh, Trlwe8{param.k, param.N}));
    std::vector<std::vector<TrlweDft24>> tmpDBNtt(param.d, std::vector<TrlweDft24>(param.dh, TrlweDft24{param.k, param.N}));
    std::vector<TrlweDft24> trlweDftRes(param.d, TrlweDft24{param.k, param.N});

    syncGadgetDecomp(tmpD, trlweInput, param);
    broadcastCRT(tmpDB, tmpD, param);

    // ntt
    for (size_t d = 0; d < param.d; d++) {
        for (size_t dh = 0; dh < param.dh; dh++) {
            for (size_t k = 0; k < param.k; k++) {
                applyNtt24(tmpDBNtt[d][dh].a[k], tmpDB[d][dh].a[k]);
            }
            applyNtt24(tmpDBNtt[d][dh].b, tmpDB[d][dh].b);
        }
    }

/*    for (size_t d = 0; d < param.d; d++) {
        auto& rgswIn = trgswDftInput[d];
        auto& rlweResA = trlweDftRes[d].a;
        auto& rlweResB = trlweDftRes[d].b;
        for (size_t l = 0; l < param.dh; l++) {
            auto& decompNttA = tmpDBNtt[d][l].a;
            auto& decompNttB = tmpDBNtt[d][l].b;
            for (size_t k = 0; k < param.k; k++) {
                auto& rgswNttA = rgswIn.trlweDftSamples[l][k].a;
                auto& rgswNttB = rgswIn.trlweDftSamples[l][k].b;
                for (size_t ka = 0; ka < param.k; ka++) {
                    for (size_t j = 0; j < param.N; j++) {
                        auto tmpAA = modMULT24(decompNttA[ka].coeffs[j], rgswNttA[ka].coeffs[j]);
                        auto tmpAB = modMULT24(decompNttA[ka].coeffs[j], rgswNttB.coeffs[j]);
                        rlweResA[ka].coeffs[j] = modADD24(rlweResA[ka].coeffs[j], tmpAA);
                        rlweResB.coeffs[j] = modADD24(rlweResB.coeffs[j], tmpAB);
                    }
                }
            }
            auto& rgswNttAk = rgswIn.trlweDftSamples[l][param.k].a;
            auto& rgswNttBk = rgswIn.trlweDftSamples[l][param.k].b;
            for (size_t ka = 0; ka < param.k; ka++) {
                for (size_t j = 0; j < param.N; j++) {
                    auto tmpBA = modMULT24(decompNttB.coeffs[j], rgswNttAk[ka].coeffs[j]);
                    rlweResA[ka].coeffs[j] = modADD24(rlweResA[ka].coeffs[j], tmpBA);
                }
            }
            for (size_t j = 0; j < param.N; j++) {
                auto tmpBB = modMULT24(decompNttB.coeffs[j], rgswNttBk.coeffs[j]);
                rlweResB.coeffs[j] = modADD24(rlweResB.coeffs[j], tmpBB);
            }
        }
    }*/

    for (size_t d = 0; d < param.d; d++) {
        auto& rgswIn = trgswDftInput[d];
        auto& rlweResA = trlweDftRes[d].a;
        auto& rlweResB = trlweDftRes[d].b;
        for (size_t l = 0; l < param.dh; l++) {
            auto& decompNttA = tmpDBNtt[d][l].a;
            auto& decompNttB = tmpDBNtt[d][l].b;
            for (size_t k = 0; k < param.k + 1; k++) {
                auto& rgswNttA = rgswIn.trlweDftSamples[l][k].a;
                auto& rgswNttB = rgswIn.trlweDftSamples[l][k].b;
                auto& decompNtt = (k < param.k) ? decompNttA[k] : decompNttB;
                for (size_t ka = 0; ka < param.k + 1; ka++) {
                    auto& rgswNtt = (ka < param.k) ? rgswNttA[ka] : rgswNttB;
                    auto& rlweRes = (ka < param.k) ? rlweResA[ka] : rlweResB;
                    for (size_t j = 0; j < param.N; j++) {
                        auto tmp = modMULT24(decompNtt.coeffs[j], rgswNtt.coeffs[j]);
                        rlweRes.coeffs[j] = modADD24(rlweRes.coeffs[j], tmp);
                    }
                }
            }
        }
    }

    // intt
    for (size_t d = 0; d < param.d; d++) {
        for (size_t k = 0; k < param.k; k++) {
            applyIntt24(output[d].a[k], trlweDftRes[d].a[k]);
        }
        applyIntt24(output[d].b, trlweDftRes[d].b);
    }
}
