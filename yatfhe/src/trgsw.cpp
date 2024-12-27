//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/crt.h"
#include "yautil/tool.h"

//todo
void trgswMPEncrypt(TrgswMP& trgswMP, const Integer mu, const YatfheParameters& param, const TrgswKey& trgswKey) {
    std::vector<TorusPolynomial> muPoly(param.k, TorusPolynomial(param.N));
    int pos = 0;
    for (size_t lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly[0].coeffs[pos] = decomposedMu;
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly[0].coeffs);
        symEncTrlweSingleSample(trgswMP.c[lvl], trgswKey.trlweKey, 0);
        muPoly[0].coeffs[pos] = mu;
        for (size_t i = 0; i < param.k; i++) {
            polynomialAddT32(trgswMP.c[lvl].a[i], trgswMP.c[lvl].a[i], muPoly[i]);
        }
    }
}

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
/*    // add the diagonal matrix (mu * G^T)_ijk to the output
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
    // ( a_0  a_1          a_k-1  b  )*/

    int pos = 0;
    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits -  (lvl + 1) * param.radixBits);
        // todo: decompose on second level
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[pos] = modAddT32(trgsw.trlweSamples[lvl][row].a[row].coeffs[pos], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[pos] = modAddT32(trgsw.trlweSamples[lvl][row].b.coeffs[pos], decomposedMu);
        }
    }
}

void trgswRotate(Trgsw& trgsw, const int rot, const YatfheParameters& param) {
    Trlwe rotT{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            rotT = trgsw.trlweSamples[lvl][row];
            trlweRotate(trgsw.trlweSamples[lvl][row], rotT, rot);
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
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = modAddT32(trgsw.trlweSamples[lvl][row].a[row].coeffs[0], decomposedMu);
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = modAddT32(trgsw.trlweSamples[lvl][row].b.coeffs[0], decomposedMu);
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
        }
    }
}

void trgswAddIntegerApproxCRT(Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = modMulT32(mu, param.w[lvl]);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = modAddT32(trgsw.trlweSamples[lvl][row].a[row].coeffs[0], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = modAddT32(trgsw.trlweSamples[lvl][row].b.coeffs[0], decomposedMu);
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

void trgswEncryptApproxCRT(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    trgswEncZero(trgsw, param, trgswKey);
    trgswAddIntegerApproxCRT(trgsw, mu, param);
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

void trgswMCRTDecomp(std::vector<Trgsw8>& out, const Trgsw& in, const YatfheParameters& param) {
    for (size_t d = 0; d < param.d; d++) {
        int64_t taoU = param.taoU[d];
        auto qd = param.qd[d];
        auto& outD = out[d];
        for (size_t l = 0; l < param.dh; l++) {
            for (size_t k1 = 0; k1 < param.k + 1; k1++) {
                auto& inA = in.trlweSamples[l][k1].a;
                auto& inB = in.trlweSamples[l][k1].b;
                auto& outA = outD.trlweSamples[l][k1].a;
                auto& outB = outD.trlweSamples[l][k1].b;
                for (size_t k2 = 0; k2 < param.k; k2++) {
                    auto& coeffOutA = outA[k2].coeffs;
                    auto& coeffInA = inA[k2].coeffs;
                    for (size_t j = 0; j < param.N; j++) {
                        coeffOutA[j] = static_cast<int8_t>(longModP(taoU * coeffInA[j], qd));
                    }
                }
                auto& coeffOutB = outB.coeffs;
                auto& coeffInB = inB.coeffs;
                for (size_t j = 0; j < param.N; j++) {
                    coeffOutB[j] = static_cast<int8_t>(longModP(taoU * coeffInB[j], qd));
                }
            }
        }
    }
}

void trgswMCRTToCRT(std::vector<Trgsw8>& trgsw, const YatfheParameters& param) {
    auto dh = param.dh;
    for (size_t d = 0; d < param.dl; d++) {
        auto ql = param.ql[d];
        int64_t taoUInv = param.taoUInv[dh + d];
        auto& inDl = trgsw[dh + d];
        for (size_t l = 0; l < param.dh; l++) {
            for (size_t k1 = 0; k1 < param.k + 1; k1++) {
                auto& inA = inDl.trlweSamples[l][k1].a;
                auto& inB = inDl.trlweSamples[l][k1].b;
                for (size_t k2 = 0; k2 < param.k; k2++) {
                    auto& coeffInA = inA[k2].coeffs;
                    for (size_t j = 0; j < param.N; j++) {
                        auto aCopy = coeffInA[j];
                        coeffInA[j] = static_cast<int8_t>(longModP(taoUInv * aCopy, ql));
                    }
                }
                auto& coeffInB = inB.coeffs;
                for (size_t j = 0; j < param.N; j++) {
                    auto bCopy = coeffInB[j];
                    coeffInB[j] = static_cast<int8_t>(longModP(taoUInv * bCopy, ql));
                }
            }
        }
    }
}

void trgswCRTRecomp(Trgsw& out, const std::vector<Trgsw8>& in, const YatfheParameters& param) {
    auto qCRT = param.qCRT;
    for (size_t l = 0; l < param.dh; l++) {
        for (size_t k1 = 0; k1 < param.k + 1; k1++) {
            auto& outA = out.trlweSamples[l][k1].a;
            for (size_t k2 = 0; k2 < param.k; k2++) {
                auto& coeffOutA = outA[k2].coeffs;
                for (size_t j = 0; j < param.N; j++) {
                    int64_t tmpA = 0l;
                    for (size_t d = 0; d < param.d; d++) {
                        auto inA = in[d].trlweSamples[l][k1].a[k2].coeffs[j];
                        tmpA += inA * param.z[d];
                    }
                    coeffOutA[j] = static_cast<Torus>(longModP(tmpA, qCRT));
                }
            }
            auto& coeffOutB = out.trlweSamples[l][k1].b.coeffs;
            for (size_t j = 0; j < param.N; j++) {
                int64_t tmpB = 0l;
                for (size_t d = 0; d < param.d; d++) {
                    auto inB = in[d].trlweSamples[l][k1].b.coeffs[j];
                    tmpB += inB * param.z[d];
                }
                coeffOutB[j] = static_cast<Torus>(longModP(tmpB, qCRT));
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
                polynomialMulAccNaiveT32(out, curr, curr2);
            }
        }
    }
}

void trgswExternalProductNtt(Trlwe& output, const TrgswDft& trgswDftInput, const Trlwe& trlweInput, const YatfheParameters& param) {
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

void trgswExternalProductApproxCRT(std::vector<Trlwe8>& output, const std::vector<Trgsw8>& trgswInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param) {
    std::vector<Trlwe8> tmpD(param.dh, Trlwe8{param.k, param.N});
    std::vector<std::vector<Trlwe8>> tmpDB(param.d, std::vector<Trlwe8>(param.dh, Trlwe8{param.k, param.N}));

    trlweApproxCRTDecomp(tmpD, trlweInput, param);
    trlweApproxCRTBroadcast(tmpDB, tmpD, param);

    for (size_t d = 0; d < param.d; d++) {
        auto& rgswIn = trgswInput[d];
        auto& rlweResA = output[d].a;
        auto& rlweResB = output[d].b;
        auto qd = param.qd[d];
        for (size_t l = 0; l < param.dh; l++) {
            auto& decompA = tmpDB[d][l].a;
            auto& decompB = tmpDB[d][l].b;
            for (size_t k = 0; k < param.k + 1; k++) {
                auto& rgswA = rgswIn.trlweSamples[l][k].a;
                auto& rgswB = rgswIn.trlweSamples[l][k].b;
                auto& decomp = (k < param.k) ? decompA[k] : decompB;
                for (size_t ka = 0; ka < param.k + 1; ka++) {
                    auto& rgsw = (ka < param.k) ? rgswA[ka] : rgswB;
                    auto& rlwe = (ka < param.k) ? rlweResA[ka] : rlweResB;
                    polynomialMulAccNaiveI8(rlwe, decomp, rgsw, qd);
                }
            }
        }
    }
}

void trgswExternalProductApproxCRTNtt(std::vector<Trlwe8>& output, const std::vector<TrgswDft24>& trgswDftInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param) {
    std::vector<Trlwe8> tmpD(param.dh, Trlwe8{param.k, param.N});
    std::vector<std::vector<Trlwe8>> tmpDB(param.d, std::vector<Trlwe8>(param.dh, Trlwe8{param.k, param.N}));
    std::vector<std::vector<TrlweDft24>> tmpDBNtt(param.d, std::vector<TrlweDft24>(param.dh, TrlweDft24{param.k, param.N}));
    std::vector<TrlweDft24> trlweDftRes(param.d, TrlweDft24{param.k, param.N});

    trlweApproxCRTDecomp(tmpD, trlweInput, param);
    trlweApproxCRTBroadcast(tmpDB, tmpD, param);

    // ntt
    for (size_t d = 0; d < param.d; d++) {
        for (size_t dh = 0; dh < param.dh; dh++) {
            for (size_t k = 0; k < param.k; k++) {
                applyNtt24(tmpDBNtt[d][dh].a[k], tmpDB[d][dh].a[k]);
            }
            applyNtt24(tmpDBNtt[d][dh].b, tmpDB[d][dh].b);
        }
    }

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
        auto qd = param.qd[d];
        for (size_t k = 0; k < param.k; k++) {
            applyIntt24(output[d].a[k], trlweDftRes[d].a[k], qd);
        }
        applyIntt24(output[d].b, trlweDftRes[d].b, qd);
    }
}

//todo
void trgswMPExternalProduct(Trlwe& output, const TrgswMP& trgswMPInput, const Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    const auto level = param.l;
    DecomposedTrlwe decomposedTrlwe{param};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    Trlwe resA{k, N};
    Trlwe resB{k, N};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& cA = trgswMPInput.c[lvl].a;
        auto& cB = trgswMPInput.c[lvl].b;
        auto& cPrimeA = trgswMPInput.cPrime[lvl].a;
        auto& cPrimeB = trgswMPInput.cPrime[lvl].b;
        auto& inA = decomposedTrlwe.rlwes[lvl].a;
        auto& inB = decomposedTrlwe.rlwes[lvl].b;
        for(size_t i = 0; i < k; i++) {
            polynomialMulAccNaiveT32(resA.a[i], inA[i], cA[i]);
            polynomialMulAccNaiveT32(resA.b, inA[i], cB);
            polynomialMulAccNaiveT32(resB.a[i], inB, cPrimeA[i]);
        }
        polynomialMulAccNaiveT32(resB.b, inB, cPrimeB);
    }
    trlweAdd(output, resB, resA);
}