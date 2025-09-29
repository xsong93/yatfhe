//
// Created by Xintong Song on 2023/12/25.
//

#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
//#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/crt.h"
#include "yatfhe/trlev.h"
#include "yautil/multi_threading.h"
#include "yautil/tool.h"

using namespace NttHexl;

void encryptTrgswMP(TrgswMP& trgswMP, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    for (size_t lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly.coeffs[pos] = decomposedMu;
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        for (size_t k = 0; k < param.k; k++) {
            symEncTrlweSingleSample(trgswMP.c[lvl][k], trgswKey.trlweKey, 0);
            for (size_t i = 0; i < param.k; i++) {
                if (i == k) {
                    addTorusPolynomial(trgswMP.c[lvl][k].a[i], trgswMP.c[lvl][k].a[i], muPoly);
                }
            }
        }
    }
}

void encryptTrgswMPMulti(TrgswMP& trgswMP, const vector<Integer>& mus, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (int j = 0; j < mus.size(); j++) {
            muPoly.coeffs[j] = mus[j] << (param.torusBits - (lvl + 1) * param.radixBits);
        }
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        for (auto k = 0; k < param.k; k++) {
            symEncTrlweSingleSample(trgswMP.c[lvl][k], trgswKey.trlweKey, 0);
            addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
        }
    }
}

void encryptTrgswMPNtt(TrgswMPDft& trgswMPDft, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    TrgswMP trgswMP{param};
    TorusPolynomial muPoly{param.N};
    for (size_t lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly.coeffs[pos] = decomposedMu;
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        applyNttForAB(trgswMPDft.cPrime[lvl], trgswMP.cPrime[lvl]);
        for (size_t k = 0; k < param.k; k++) {
            symEncTrlweSingleSample(trgswMP.c[lvl][k], trgswKey.trlweKey, 0);
            addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
            applyNttForAB(trgswMPDft.c[lvl][k], trgswMP.c[lvl][k]);
        }
    }
}

void encryptTrgswMPFixedNoiseNtt(TrgswMPDft& trgswMPDft, const Integer mu, const TrgswKey& trgswKey, const int pos, const Torus noise, const YatfheParameters& param) {
    TrgswMP trgswMP{param};
    TorusPolynomial muPoly{param.N};
    for (size_t lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly.coeffs[pos] = decomposedMu;
        symEncTrlweMultiSampleFixedNoise(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs, noise);
        applyNttForAB(trgswMPDft.cPrime[lvl], trgswMP.cPrime[lvl]);
        for (size_t k = 0; k < param.k; k++) {
            symEncTrlweSingleSampleFixedNoise(trgswMP.c[lvl][k], trgswKey.trlweKey, 0, noise);
            addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
            applyNttForAB(trgswMPDft.c[lvl][k], trgswMP.c[lvl][k]);
        }
    }
}

void encryptTrgswMPMultiNtt(TrgswMPDft& trgswMPDft, const vector<Integer>& mus, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TrgswMP trgswMP{param};
    TorusPolynomial muPoly{param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (int j = 0; j < mus.size(); j++) {
            muPoly.coeffs[j] = mus[j] << (param.torusBits - (lvl + 1) * param.radixBits);
        }
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        applyNttForAB(trgswMPDft.cPrime[lvl], trgswMP.cPrime[lvl]);
        for (auto k = 0; k < param.k; k++) {
            symEncTrlweSingleSample(trgswMP.c[lvl][k], trgswKey.trlweKey, 0);
            addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
            applyNttForAB(trgswMPDft.c[lvl][k], trgswMP.c[lvl][k]);
        }
    }
}

void encryptLowTrgswMP(TrgswMP& trgswMP, const Integer mu, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    int pos = 0;
    muPoly.coeffs[pos] = mu;
    symEncTrlweMultiSample(trgswMP.cPrime[0], trgswKey.trlweKey, muPoly.coeffs);
    for (size_t k = 0; k < param.k; k++) {
        symEncTrlweSingleSample(trgswMP.c[0][k], trgswKey.trlweKey, 0);
        for (size_t i = 0; i < param.k; i++) {
            if (i == k) {
                addTorusPolynomial(trgswMP.c[0][k].a[i], trgswMP.c[0][k].a[i], muPoly);
            }
        }
    }
}

void encryptLowTrgswMPNtt(TrgswMP& trgswMP, TrgswMPDft& trgswMPDft, const Integer mu, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    int pos = 0;
    muPoly.coeffs[pos] = mu;
    symEncTrlweMultiSample(trgswMP.cPrime[0], trgswKey.trlweKey, muPoly.coeffs);
    applyNttForAB(trgswMPDft.cPrime[0], trgswMP.cPrime[0]);
    for (size_t k = 0; k < param.k; k++) {
        symEncTrlweSingleSample(trgswMP.c[0][k], trgswKey.trlweKey, 0);
        for (size_t i = 0; i < param.k; i++) {
            if (i == k) {
                addTorusPolynomial(trgswMP.c[0][k].a[i], trgswMP.c[0][k].a[i], muPoly);
            }
        }
        applyNttForAB(trgswMPDft.c[0][k], trgswMP.c[0][k]);
    }
}

// trgsw(0): [trlwe(0)]  (k+1)l rows
void encZeroTrgsw(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSample(trgsw.trlweSamples[lvl][row], trgswKey.trlweKey, 0);
        }
    }
}

// output += mu * G^T
void addIntegerToTrgsw(Trgsw& trgsw, const Integer mu, const int pos, const YatfheParameters& param) {
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

    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits -  (lvl + 1) * param.radixBits);
        // todo: decompose on second level
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[pos] = addTorus(trgsw.trlweSamples[lvl][row].a[row].coeffs[pos], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[pos] = addTorus(trgsw.trlweSamples[lvl][row].b.coeffs[pos], decomposedMu);
        }
    }
}

void rotateTrgsw(Trgsw& trgsw, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    Trlwe rotT{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            rotT = trgsw.trlweSamples[lvl][row];
            rotateTrlwe(trgsw.trlweSamples[lvl][row], rotT, rot);
        }
    }
}

void rotateTrgswNtt(TrgswDft& trgswDft, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    TrlweDft rotT{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            rotT = trgswDft.trlweDftSamples[lvl][row];
            rotateTrlweNtt(trgswDft.trlweDftSamples[lvl][row], rotT, rot);
        }
    }
}

void rotateTrgswMP(TrgswMP& trgswMP, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    Trlwe rotT1{param.k, param.N};
    Trlwe rotT2{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        rotT1 = trgswMP.cPrime[lvl];
        rotateTrlwe(trgswMP.cPrime[lvl], rotT1, rot);
        for (auto row = 0; row < param.k; row++) {
            rotT2 = trgswMP.c[lvl][row];
            rotateTrlwe(trgswMP.c[lvl][row], rotT2, rot);
        }
    }
}

void rotateTrgswMPNtt(TrgswMPDft& trgswMP, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    TrlweDft rotT1{param.k, param.N};
    TrlweDft rotT2{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        rotT1 = trgswMP.cPrime[lvl];
        rotateTrlweNtt(trgswMP.cPrime[lvl], rotT1, rot);
        for (auto row = 0; row < param.k; row++) {
            rotT2 = trgswMP.c[lvl][row];
            rotateTrlweNtt(trgswMP.c[lvl][row], rotT2, rot);
        }
    }
}

// trgsw(0): [trlwe(0)]  (k+1)l rows
void encZeroTrgswNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSampleNtt(trgsw.trlweSamples[lvl][row], trgswDft.trlweDftSamples[lvl][row], trgswKey.trlweKey, 0);
        }
    }
}

// output += mu * G^T
void addIntegerToTrgswNtt(TrgswDft& trgswDft, Trgsw& trgsw, const Integer mu, const int pos, const YatfheParameters& param) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = mu << (param.torusBits -  (lvl + 1) * param.radixBits);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = addTorus(trgsw.trlweSamples[lvl][row].a[row].coeffs[pos], decomposedMu);
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = addTorus(trgsw.trlweSamples[lvl][row].b.coeffs[pos], decomposedMu);
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
        }
    }
}

void trgswAddIntegerApproxCRT(Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
    for (auto lvl = 0; lvl < param.l; lvl++) {
        auto decomposedMu = multTorus(mu, param.w[lvl]);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = addTorus(trgsw.trlweSamples[lvl][row].a[row].coeffs[0], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = addTorus(trgsw.trlweSamples[lvl][row].b.coeffs[0], decomposedMu);
        }
    }
}

void encryptTrgsw(Trgsw& trgsw, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    encZeroTrgsw(trgsw, param, trgswKey);
    addIntegerToTrgsw(trgsw, mu, pos, param);
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
void encryptTrgswNtt(Trgsw& trgsw, TrgswDft& trgswDft, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    encZeroTrgswNtt(trgsw, trgswDft, param, trgswKey);
    addIntegerToTrgswNtt(trgswDft, trgsw, mu, pos, param);
}

void encryptTrgswApproxCRT(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    encZeroTrgsw(trgsw, param, trgswKey);
    trgswAddIntegerApproxCRT(trgsw, mu, param);
}

// To decrypt, it is sufficient to decrypt the last GLev ciphertext, which is a GLev encryption of m/B^l.
Integer decryptTrgsw(const Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRounding(tmp, trgsw.trlweSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweStdDev, param.torusBits - param.radixBits);
}

Integer decryptTrgswNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRoundingNtt(tmp, trgswDft.trlweDftSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweStdDev, param.torusBits - param.radixBits);
}

void decryptTrgswMP(IntPolynomial& res, const TrgswMP& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey, const bool isDecC) {
    const auto firstLevel = 0;
    TorusPolynomial tmp{param.N};
    if (!isDecC) {
        symDecTrlweWoRounding(tmp, trgsw.cPrime[firstLevel], trgswKey.trlweKey);
    } else {
        symDecTrlweWoRounding(tmp, trgsw.c[firstLevel][0], trgswKey.trlweKey);
    }
    for (auto i = 0; i < param.N; i++) {
        res.coeffs[i] = roundErrorForShiftedTorus(tmp.coeffs[i], param.rlweStdDev, param.torusBits - param.radixBits);
    }
}

void decryptTrgswMPNtt(IntPolynomial& res, const TrgswMPDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey, const bool isDecC) {
    const auto firstLevel = 0;
    TorusPolynomial tmp{param.N};
    if (!isDecC) {
        symDecTrlweWoRoundingNtt(tmp, trgswDft.cPrime[firstLevel], trgswKey.trlweKey);
    } else {
        symDecTrlweWoRoundingNtt(tmp, trgswDft.c[firstLevel][0], trgswKey.trlweKey);
    }
    for (auto i = 0; i < param.N; i++) {
        res.coeffs[i] = roundErrorForShiftedTorus(tmp.coeffs[i], param.rlweStdDev, param.torusBits - param.radixBits);
    }
}

void decompTrgswMcrt(std::vector<Trgsw8>& out, const Trgsw& in, const YatfheParameters& param) {
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

void trgswMcrtToCrt(std::vector<Trgsw8>& trgsw, const YatfheParameters& param) {
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

void recompTrgswCrt(Trgsw& out, const std::vector<Trgsw8>& in, const YatfheParameters& param) {
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

void externalProductTrgsw(Trlwe& output, const Trgsw& trgswInput, const Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level = trgswInput.l;
    DecomposedTrlwe decomposedTrlwe {param};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    // <Decomp(B), C_k^bar> + Σ_(0,k-1)<Decomp(A_i), C_i^bar>
    // https://www.zama.ai/post/tfhe-deep-dive-part-3
    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto& curr = (col < k) ? decomposedTrlwe.trlwes[lvl].a[col] : decomposedTrlwe.trlwes[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                auto& out = (col2 < k) ? output.a[col2] : output.b;
                auto& curr2 = (col2 < k) ? trgswInput.trlweSamples[lvl][col].a[col2]
                                         : trgswInput.trlweSamples[lvl][col].b;
                multTorusPolynomialAcc(out, curr, curr2);
            }
        }
    }
}

void externalProductTrgswNtt(Trlwe& output, const TrgswDft& trgswDftInput, const Trlwe& trlweInput, const int level, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto N = trlweInput.b.N;
    TrlweDft trlweDftRes {k, N};
    DecomposedTrlwe decomposedTrlwe {param};
    DecomposedTrlweDft decomposedTrlweDft {param, level};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);  // 8 * 2

//#pragma omp parallel for
    for (auto i = 0; i < level; i++) {
        applyNttForAB(decomposedTrlweDft.rlweDfts[i], decomposedTrlwe.trlwes[i]);
    }

//#pragma omp parallel for collapse(2) private(out)
    for (auto lvl = 0; lvl < level; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto& curr = (col < k) ? decomposedTrlweDft.rlweDfts[lvl].a[col] : decomposedTrlweDft.rlweDfts[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                auto& out = (col2 < k) ? trlweDftRes.a[col2] : trlweDftRes.b;
                auto& curr2 = (col2 < k) ? trgswDftInput.trlweDftSamples[lvl][col].a[col2]
                                         : trgswDftInput.trlweDftSamples[lvl][col].b;
                NttPolynomial tmp{N};
//                EltwiseMultMod(tmp.coeffs.data(), curr.coeffs.data(), curr2.coeffs.data(), N, param.qNtt, 1);
//                EltwiseAddMod(out.coeffs.data(), out.coeffs.data(), tmp.coeffs.data(), N, param.qNtt);
//                calModularInnerProductNttHexl(out, curr, curr2, N, param.qNtt);
                calModularInnerProductNtt(out, curr, curr2);
            }
        }
    }
    applyInttForAB(output, trlweDftRes);
}

void externalProductTrgswApproxCrt(std::vector<Trlwe8>& output, const std::vector<Trgsw8>& trgswInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param) {
    std::vector<Trlwe8> tmpD(param.dh, Trlwe8{param.k, param.N});
    std::vector<std::vector<Trlwe8>> tmpDB(param.d, std::vector<Trlwe8>(param.dh, Trlwe8{param.k, param.N}));

    decompTrlweApproxCrt(tmpD, trlweInput, param);
    broadcastTrlweApproxCrt(tmpDB, tmpD, param);

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
                    multInt8PolynomialAcc(rlwe, decomp, rgsw, qd);
                }
            }
        }
    }
}

void externalProductTrgswApproxCrtNtt(std::vector<Trlwe8>& output, const std::vector<TrgswDft24>& trgswDftInput, const std::vector<Trlwe8>& trlweInput, const YatfheParameters& param) {
    std::vector<Trlwe8> tmpD(param.dh, Trlwe8{param.k, param.N});
    std::vector<std::vector<Trlwe8>> tmpDB(param.d, std::vector<Trlwe8>(param.dh, Trlwe8{param.k, param.N}));
    std::vector<std::vector<TrlweDft24>> tmpDBNtt(param.d, std::vector<TrlweDft24>(param.dh, TrlweDft24{param.k, param.N}));
    std::vector<TrlweDft24> trlweDftRes(param.d, TrlweDft24{param.k, param.N});

    decompTrlweApproxCrt(tmpD, trlweInput, param);
    broadcastTrlweApproxCrt(tmpDB, tmpD, param);

    // ntt
    for (size_t d = 0; d < param.d; d++) {
        for (size_t dh = 0; dh < param.dh; dh++) {
            for (size_t k = 0; k < param.k; k++) {
                NttNative24::applyNtt(tmpDBNtt[d][dh].a[k], tmpDB[d][dh].a[k]);
            }
            NttNative24::applyNtt(tmpDBNtt[d][dh].b, tmpDB[d][dh].b);
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
                        auto tmp = NttNative24::modMULT(decompNtt.coeffs[j], rgswNtt.coeffs[j]);
                        rlweRes.coeffs[j] = NttNative24::modADD(rlweRes.coeffs[j], tmp);
                    }
                }
            }
        }
    }

    // intt
    for (size_t d = 0; d < param.d; d++) {
        auto qd = param.qd[d];
        for (size_t k = 0; k < param.k; k++) {
            NttNative24::applyIntt(output[d].a[k], trlweDftRes[d].a[k], qd);
        }
        NttNative24::applyIntt(output[d].b, trlweDftRes[d].b, qd);
    }
}

void externalProductTrgswMP(Trlwe& output, const TrgswMP& trgswMPInput, const Trlwe& trlweInput, const int level, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    DecomposedTrlwe decomposedTrlwe{param};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    Trlwe resA{k, N};
    Trlwe resB{k, N};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[lvl];
        auto& cPrimeA = trgswMPInput.cPrime[lvl].a;
        auto& cPrimeB = trgswMPInput.cPrime[lvl].b;
        auto& inA = decomposedTrlwe.trlwes[lvl].a;
        auto& inB = decomposedTrlwe.trlwes[lvl].b;
        for(size_t i = 0; i < k; i++) {
            auto& ciA = c[i].a;
            auto& ciB = c[i].b;
            for (size_t i2 = 0; i2 < k; i2++) {
                multTorusPolynomialAcc(resA.a[i2], inA[i], ciA[i2]);
            }
            multTorusPolynomialAcc(resA.b, inA[i], ciB);
            multTorusPolynomialAcc(resB.a[i], inB, cPrimeA[i]);
        }
        multTorusPolynomialAcc(resB.b, inB, cPrimeB);
    }
    addTrlwe(output, resB, resA);
}

void externalProductTrgswMPNtt(Trlwe& output, const TrgswMPDft& trgswMPInput, const Trlwe& trlweInput, const int level, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    DecomposedTrlwe decomposedTrlwe{param};
    DecomposedTrlweDft decomposedTrlweDft{param, level};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    for (auto i = 0; i < level; i++) {
        applyNttForAB(decomposedTrlweDft.rlweDfts[i], decomposedTrlwe.trlwes[i]);
    }

    TrlweDft resA{k, N};
    TrlweDft resB{k, N};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[lvl];
        auto& cPrimeA = trgswMPInput.cPrime[lvl].a;
        auto& cPrimeB = trgswMPInput.cPrime[lvl].b;
        auto& inA = decomposedTrlweDft.rlweDfts[lvl].a;
        auto& inB = decomposedTrlweDft.rlweDfts[lvl].b;
        for(size_t i = 0; i < k; i++) {
            auto& ciA = c[i].a;
            auto& ciB = c[i].b;
            for (size_t i2 = 0; i2 < k; i2++) {
                calModularInnerProductNtt(resA.a[i2], inA[i], ciA[i2]);
            }
            calModularInnerProductNtt(resA.b, inA[i], ciB);
            calModularInnerProductNtt(resB.a[i], inB, cPrimeA[i]);
        }
        calModularInnerProductNtt(resB.b, inB, cPrimeB);
    }
    TrlweDft tmp{k, N};
    addTrlweNtt(tmp, resB, resA);
    applyInttForAB(output, tmp);
}

void externalProductTrgswMPNttInPlace(Trlwe& acc, const TrgswMPDft& trgswMPInput, const int level, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    DecomposedTrlwe decomposedTrlwe{param};
    DecomposedTrlweDft decomposedTrlweDft{param, level};
    gadgetDecomposeTrlwe(decomposedTrlwe, acc, param);

    for (auto i = 0; i < level; i++) {
        applyNttForAB(decomposedTrlweDft.rlweDfts[i], decomposedTrlwe.trlwes[i]);
    }

    TrlweDft resA{k, N};
    TrlweDft resB{k, N};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[lvl];
        auto& cPrimeA = trgswMPInput.cPrime[lvl].a;
        auto& cPrimeB = trgswMPInput.cPrime[lvl].b;
        auto& inA = decomposedTrlweDft.rlweDfts[lvl].a;
        auto& inB = decomposedTrlweDft.rlweDfts[lvl].b;
        for(size_t i = 0; i < k; i++) {
            auto& ciA = c[i].a;
            auto& ciB = c[i].b;
            for (size_t i2 = 0; i2 < k; i2++) {
                calModularInnerProductNtt(resA.a[i2], inA[i], ciA[i2]);
            }
            calModularInnerProductNtt(resA.b, inA[i], ciB);
            calModularInnerProductNtt(resB.a[i], inB, cPrimeA[i]);
        }
        calModularInnerProductNtt(resB.b, inB, cPrimeB);
    }
    TrlweDft tmp{k, N};
    addTrlweNtt(tmp, resB, resA);
    applyInttForAB(acc, tmp);
}


void externalProductTrgswMPNttMT(Trlwe& output, const TrgswMPDft& trgswMPInput, const Trlwe& trlweInput, const int level, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    DecomposedTrlwe decomposedTrlwe{param};
    DecomposedTrlweDft decomposedTrlweDft{param, level};
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(level);
    vector resAV(level, TrlweDft{k, N});
    vector resBV(level, TrlweDft{k, N});

    for (auto lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[lvl];
        auto& cPrimeA = trgswMPInput.cPrime[lvl].a;
        auto& cPrimeB = trgswMPInput.cPrime[lvl].b;
        auto& in = decomposedTrlwe.trlwes[lvl];
        auto& resA = resAV[lvl];
        auto& resB = resBV[lvl];
        futures.emplace_back(pool.enqueue([k, N, &c, &cPrimeA, &cPrimeB, &in, &resA, &resB] {
            TrlweDft inDft{k, N};
            applyNttForAB(inDft, in);
            auto& inA = inDft.a;
            auto& inB = inDft.b;
            for(size_t i = 0; i < k; i++) {
                auto& ciA = c[i].a;
                auto& ciB = c[i].b;
                for (size_t i2 = 0; i2 < k; i2++) {
                    calModularInnerProductNtt(resA.a[i2], inA[i], ciA[i2]);
                }
                calModularInnerProductNtt(resA.b, inA[i], ciB);
                calModularInnerProductNtt(resB.a[i], inB, cPrimeA[i]);
            }
            calModularInnerProductNtt(resB.b, inB, cPrimeB);
        }));
    }
    for (auto& f : futures) {
        f.wait();
    }
    TrlweDft tmp{k, N};
    for (size_t lvl = 0; lvl < level; lvl++) {
        addTrlweNtt(tmp, tmp, resAV[lvl]);
        addTrlweNtt(tmp, tmp, resBV[lvl]);
    }
    applyInttForAB(output, tmp);
}

void externalProductTrgswMPDecomp(DecomposedTrlwe& output, const TrgswMP& trgswMPInput, const DecomposedTrlwe& trlweInput, const YatfheParameters& param) {
    const auto k = param.k;
    const auto level = param.l;

    DecomposedTrlwe resA{param};
    DecomposedTrlwe resB{param};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[0];
        auto& cPrimeA = trgswMPInput.cPrime[0].a;
        auto& cPrimeB = trgswMPInput.cPrime[0].b;
        auto& inA = trlweInput.trlwes[lvl].a;
        auto& inB = trlweInput.trlwes[lvl].b;
        for(size_t i = 0; i < k; i++) {
            auto& ciA = c[i].a;
            auto& ciB = c[i].b;
            for (size_t i2 = 0; i2 < k; i2++) {
                multTorusPolynomialAcc(resA.trlwes[lvl].a[i2], inA[i], ciA[i2]);
            }
            multTorusPolynomialAcc(resA.trlwes[lvl].b, inA[i], ciB);
            multTorusPolynomialAcc(resB.trlwes[lvl].a[i], inB, cPrimeA[i]);
        }
        multTorusPolynomialAcc(resB.trlwes[lvl].b, inB, cPrimeB);
    }
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& out = output.trlwes[lvl];
        auto& a = resA.trlwes[lvl];
        auto& b = resB.trlwes[lvl];
        for (auto i = 0; i < param.k; i++) {
            addTorusPolynomial(out.a[i], a.a[i], b.a[i]);
        }
        addTorusPolynomial(out.b, a.b, b.b);
    }
}

void externalProductTrgswMPDecompNtt(DecomposedTrlweDft& output, const TrgswMPDft& trgswMPInput, const DecomposedTrlweDft& trlweInput, const YatfheParameters& param) {
    const auto k = param.k;
    const auto level = param.l;

    DecomposedTrlweDft resA{param, level};
    DecomposedTrlweDft resB{param, level};
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& c = trgswMPInput.c[0];
        auto& cPrimeA = trgswMPInput.cPrime[0].a;
        auto& cPrimeB = trgswMPInput.cPrime[0].b;
        auto& inA = trlweInput.rlweDfts[lvl].a;
        auto& inB = trlweInput.rlweDfts[lvl].b;
        for(size_t i = 0; i < k; i++) {
            auto& ciA = c[i].a;
            auto& ciB = c[i].b;
            for (size_t i2 = 0; i2 < k; i2++) {
                calModularInnerProductNtt(resA.rlweDfts[lvl].a[i2], inA[i], ciA[i2]);
            }
            calModularInnerProductNtt(resA.rlweDfts[lvl].b, inA[i], ciB);
            calModularInnerProductNtt(resB.rlweDfts[lvl].a[i], inB, cPrimeA[i]);
        }
        calModularInnerProductNtt(resB.rlweDfts[lvl].b, inB, cPrimeB);
    }
    for (size_t lvl = 0; lvl < level; lvl++) {
        auto& out = output.rlweDfts[lvl];
        auto& a = resA.rlweDfts[lvl];
        auto& b = resB.rlweDfts[lvl];
        addTrlweNtt(out, a, b);
    }
}

void generalExternalProductTrgswMPNtt(Trlev& output, const TrgswMPDft& input1, const Trlev& input2, const int level, const YatfheParameters& param) {
    const auto L = level;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(L);
    for (auto l = 0; l < L; l++) {
        futures.emplace_back(pool.enqueue([&output, &input1, &input2, &param, l] {
            externalProductTrgswMPNtt(output.trlwes[l], input1, input2.trlwes[l], param.l, param);
        }));
    }
    for (auto& f : futures) {
        f.wait();
    }
}

void internalProductTrgswMP(TrgswMP& output, const TrgswMP& input1, const TrgswMP& input2, const int level, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = level;
    for (size_t l = 0; l < L; l++) {
        externalProductTrgswMP(output.cPrime[l], input1, input2.cPrime[l], param.l, param);
        for (size_t k = 0; k < K; k++) {
            externalProductTrgswMP(output.c[l][k], input1, input2.c[l][k], param.l, param);
        }
    }
}

void internalProductTrgswMPNtt(TrgswMP& output, const TrgswMP& input1, const TrgswMPDft& input2, const int level, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = level;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(L + K * L);
    for (size_t l = 0; l < L; l++) {
        futures.emplace_back(pool.enqueue([&output, &input1, &input2, &param, l] {
            externalProductTrgswMPNtt(output.cPrime[l], input2, input1.cPrime[l], param.l, param);
        }));
    }
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K; k++) {
            futures.emplace_back(pool.enqueue([&output, &input1, &input2, &param, l, k] {
                 externalProductTrgswMPNtt(output.c[l][k], input2, input1.c[l][k], param.l, param);
            }));
        }
    }
    for (auto& f : futures) {
        f.wait();
    }
}

void internalProductTrgswMPNtt(TrgswMPDft& output, const TrgswMP& input1, const TrgswMPDft& input2, const int level, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = level;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(L + K * L);
    for (size_t l = 0; l < L; l++) {
        auto& cPrimeL = output.cPrime[l];
        futures.emplace_back(pool.enqueue([&cPrimeL, &input1, &input2, &param, l] {
            Trlwe tmp{param.k, param.N};
            externalProductTrgswMPNtt(tmp, input2, input1.cPrime[l], param.l, param);
            applyNttForAB(cPrimeL, tmp);
        }));
    }
    for (size_t l = 0; l < L; l++) {
        for (size_t k = 0; k < K; k++) {
            auto& c = output.c[l][k];
            futures.emplace_back(pool.enqueue([&c, &input1, &input2, &param, l, k] {
                Trlwe tmp{param.k, param.N};
                externalProductTrgswMPNtt(tmp, input2, input1.c[l][k], param.l, param);
                applyNttForAB(c, tmp);
            }));
        }
    }
    for (auto& f : futures) {
        f.wait();
    }

// #pragma omp parallel for simd collapse(2) schedule(guided)
//     for (size_t l = 0; l < L; l++) {
//         for (size_t k = 0; k < K; k++) {
//             if (k == 0) {
//                 externalProductTrgswMPNtt(output.cPrime[l], input2, input1.cPrime[l], param);
//             }
//             externalProductTrgswMPNtt(output.c[l][k], input2, input1.c[l][k], param);
//         }
//     }
}

void internalProductAsymTrgswMPNtt(TrgswMPDft& output, const TrgswMPDft& input1, const Trlev& input2, const TrlevDft& sSquare, const int level, const YatfheParameters& param) {
    const auto L = level;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(L);
    for (auto l = 0; l < L; l++) {
        auto& cPrimeL = output.cPrime[l];
        auto& cL = output.c[l];
        auto& in2L = input2.trlwes[l];
        futures.emplace_back(pool.enqueue([&cPrimeL, &cL, &input1, &in2L, &sSquare, &param] {
            Trlwe tmp{param.k, param.N};
            vector tmpC(param.k, Trlwe{param.k, param.N});
            externalProductTrgswMPNtt(tmp, input1, in2L, param.l, param);
            switchTrlevToTrgswNtt(tmpC, tmp, sSquare, param);
            applyNttForAB(cPrimeL, tmp);
            for (auto i = 0; i < param.k; i++) {
                applyNttForAB(cL[i], tmpC[i]);
            }
        }));
    }
    for (auto& f : futures) {
        f.wait();
    }
}

void switchTrlevToTrgswNtt(vector<Trlwe>& c, const Trlwe& cPrime, const TrlevDft& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l;
    const auto N = param.N;
    auto& cPrimeA = cPrime.a;
    auto& cPrimeB = cPrime.b;

    vector decomp(L, Trlwe{K, N});
    for (auto row = 0; row < K; row++) {
        auto& currIn = cPrimeA[row];
        for (auto j = 0; j < N; j++) {
            DecomposedData d {L};
            gadgetDecompose(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < L; lvl++) {
                auto& currOut = decomp[lvl].a[row];
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }

    vector cDft(K, TrlweDft{K, N});
    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlweDfts[l];
        auto& decompL = decomp[l];
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = cDft[k1].a;
            auto& cB = cDft[k1].b;
            auto& sA = s2.a;
            NttPolynomial aDft{N};
            applyNtt(aDft, decompL.a[k1]);
            for (auto k2 = 0; k2 < K; k2++) {
                calModularInnerProductNtt(cA[k2], aDft, sA[k2]);
            }
            calModularInnerProductNtt(cB, aDft, s2.b);
        }
    }

    for (auto k1 = 0; k1 < K; k1++) {
        applyInttForAB(c[k1], cDft[k1]);
        for (auto k2 = 0; k2 < K; k2++) {
            addTorusPolynomial(c[k1].a[k2], c[k1].a[k2], cPrimeB);
        }
    }
}

void switchTrlevToTrgsw(vector<Trlwe>& c, const Trlwe& cPrime, const Trlev& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l;
    const auto N = param.N;
    auto& cPrimeA = cPrime.a;
    auto& cPrimeB = cPrime.b;

    vector decomp(L, Trlwe{K, N});

    for (auto row = 0; row < K; row++) {
        auto& currIn = cPrimeA[row];
        for (auto j = 0; j < N; j++) {
            DecomposedData d {L};
            gadgetDecompose(d, currIn.coeffs[j], param);
            for (auto lvl = 0; lvl < L; lvl++) {
                auto& currOut = decomp[lvl].a[row];
                currOut.coeffs[j] = d.value[lvl] * d.sign;
            }
        }
    }

    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlwes[l];
        auto& decompL = decomp[l];
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = c[k1].a;
            auto& cB = c[k1].b;
            auto& sA = s2.a;
            for (auto k2 = 0; k2 < K; k2++) {
                multTorusPolynomialAcc(cA[k2], decompL.a[k1], sA[k2]);
            }
            multTorusPolynomialAcc(cB, decompL.a[k1], s2.b);
        }
    }

    for (auto k1 = 0; k1 < K; k1++) {
        for (auto k2 = 0; k2 < K; k2++) {
            addTorusPolynomial(c[k1].a[k2], c[k1].a[k2], cPrimeB);
        }
    }
}