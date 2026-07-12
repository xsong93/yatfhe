//
// Created by Xintong Song on 2023/12/25.
//

#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/crt.h"
#include "yatfhe/trlev.h"
#include "yautil/multi_threading.h"

using namespace NttHexl;

void encryptTrgswMP(TrgswMP& trgswMP, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    TrgswMPDft trgswMPDft{param, trgswMP.l, trgswMP.isHalf};
    TorusPolynomial muPoly{param.N};
    NttPolynomial myPolyNtt{param.N};
    for (size_t lvl = 0; lvl < trgswMPDft.l; lvl++) {
        auto decomposedMu = static_cast<Torus>(mu) << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly.coeffs[pos] = decomposedMu;
        applyNtt(myPolyNtt, muPoly);
        symEncTrlweMultiSampleNtt(trgswMP.cPrime[lvl], trgswMPDft.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        if (!trgswMPDft.isHalf) {
            for (size_t k = 0; k < param.k; k++) {
                symEncTrlweSingleSampleNtt(trgswMP.c[lvl][k], trgswMPDft.c[lvl][k], trgswKey.trlweKey, 0, 0);
                addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
            }
        }
    }
}

void encryptTrgswMPNtt(TrgswMPDft& trgswMPDft, const Integer mu, const TrgswKey& trgswKey, const int pos, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    NttPolynomial myPolyNtt{param.N};
    for (size_t lvl = 0; lvl < trgswMPDft.l; lvl++) {
        auto decomposedMu = static_cast<Torus>(mu) << (param.torusBits - (lvl + 1) * param.radixBits);
        muPoly.coeffs[pos] = decomposedMu;
        applyNtt(myPolyNtt, muPoly);
        symEncTrlweMultiSampleSimple(trgswMPDft.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        if (!trgswMPDft.isHalf) {
            for (size_t k = 0; k < param.k; k++) {
                symEncTrlweSingleSampleNttSimple(trgswMPDft.c[lvl][k], trgswKey.trlweKey, 0, 0);
                addNttPolynomial(trgswMPDft.c[lvl][k].a[k], trgswMPDft.c[lvl][k].a[k], myPolyNtt);
            }
        }
    }
}

void encryptTrgswMPMulti(TrgswMP& trgswMP, const vector<Integer>& mus, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    for (auto lvl = 0; lvl < trgswMP.l; lvl++) {
        for (int j = 0; j < mus.size(); j++) {
            muPoly.coeffs[j] = static_cast<Torus>(mus[j]) << (param.torusBits - (lvl + 1) * param.radixBits);
        }
        symEncTrlweMultiSample(trgswMP.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        if (!trgswMP.isHalf) {
            for (auto k = 0; k < param.k; k++) {
                symEncTrlweSingleSample(trgswMP.c[lvl][k], trgswKey.trlweKey, 0, 0);
                addTorusPolynomial(trgswMP.c[lvl][k].a[k], trgswMP.c[lvl][k].a[k], muPoly);
            }
        }
    }
}

void encryptTrgswMPMultiNtt(TrgswMPDft& trgswMPDft, const vector<Integer>& mus, const TrgswKey& trgswKey, const YatfheParameters& param) {
    TorusPolynomial muPoly{param.N};
    NttPolynomial myPolyNtt{param.N};
    for (auto lvl = 0; lvl < trgswMPDft.l; lvl++) {
        for (int j = 0; j < mus.size(); j++) {
            muPoly.coeffs[j] = static_cast<Torus>(mus[j]) << (param.torusBits - (lvl + 1) * param.radixBits);
        }
        applyNtt(myPolyNtt, muPoly);
        symEncTrlweMultiSampleSimple(trgswMPDft.cPrime[lvl], trgswKey.trlweKey, muPoly.coeffs);
        if (!trgswMPDft.isHalf) {
            for (auto k = 0; k < param.k; k++) {
                symEncTrlweSingleSampleNttSimple(trgswMPDft.c[lvl][k], trgswKey.trlweKey, 0, 0);
                addNttPolynomial(trgswMPDft.c[lvl][k].a[k], trgswMPDft.c[lvl][k].a[k], myPolyNtt);
            }
        }
    }
}

void rotateTrgsw(Trgsw& trgsw, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    Trlwe rotT{param.k, param.N};
    for (auto lvl = 0; lvl < trgsw.l; lvl++) {
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
    for (auto lvl = 0; lvl < trgswDft.l; lvl++) {
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
    for (auto lvl = 0; lvl < trgswMP.l; lvl++) {
        rotT1 = trgswMP.cPrime[lvl];
        rotateTrlwe(trgswMP.cPrime[lvl], rotT1, rot);
        for (auto row = 0; row < param.k; row++) {
            rotT2 = trgswMP.c[lvl][row];
            rotateTrlwe(trgswMP.c[lvl][row], rotT2, rot);
        }
    }
}

void rotateTrgswMPNtt(TrgswMPDft& out, const TrgswMPDft& in, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    for (auto lvl = 0; lvl < out.l; lvl++) {
        rotateTrlweNtt(out.cPrime[lvl], in.cPrime[lvl], rot);
        for (auto row = 0; row < param.k; row++) {
            rotateTrlweNtt(out.c[lvl][row], in.c[lvl][row], rot);
        }
    }
}

void rotateTrgswMPMinusOneNtt(TrgswMPDft& out, const TrgswMPDft& in, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    for (auto lvl = 0; lvl < out.l; lvl++) {
        rotateTrlweMinusOneNtt(out.cPrime[lvl], in.cPrime[lvl], rot);
        for (auto row = 0; row < param.k; row++) {
            rotateTrlweMinusOneNtt(out.c[lvl][row], in.c[lvl][row], rot);
        }
    }
}

void rotateTrgswMPMinusOneBPlusOneNtt(TrgswMPDft& out, const TrgswMPDft& in, const int rot, const YatfheParameters& param) {
  if (rot % (param.N * 2) == 0) {
      return;
  }
  const auto q = NttHexl::getNttHexl().GetModulus();
  for (auto lvl = 0; lvl < out.l; lvl++) {
      const auto g_l = static_cast<uint64_t>(
          static_cast<Torus>(1) << (param.torusBits - (lvl + 1) * param.radixBits));
      rotateTrlweMinusOneNtt(out.cPrime[lvl], in.cPrime[lvl], rot);
      for (auto& coeff : out.cPrime[lvl].b.coeffs) {
          coeff += g_l;
          if (coeff >= q) coeff -= q;
      }
      for (auto row = 0; row < param.k; row++) {
          rotateTrlweMinusOneNtt(out.c[lvl][row], in.c[lvl][row], rot);
          for (auto& coeff : out.c[lvl][row].a[row].coeffs) {
              coeff += g_l;
              if (coeff >= q) coeff -= q;
          }
      }
  }
}

// trgsw(0): [trlwe(0)]  (k+1)l rows
void encZeroTrgsw(Trgsw& trgsw, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < trgsw.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSample(trgsw.trlweSamples[lvl][row], trgswKey.trlweKey, 0, 0);
        }
    }
}

void encZeroTrgswNtt(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    for (auto lvl = 0; lvl < trgsw.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            symEncTrlweSingleSampleNtt(trgsw.trlweSamples[lvl][row], trgswDft.trlweDftSamples[lvl][row], trgswKey.trlweKey, 0, 0);
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

    for (auto lvl = 0; lvl < trgsw.l; lvl++) {
        auto decomposedMu = static_cast<Torus>(mu) << (param.torusBits -  (lvl + 1) * param.radixBits);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[pos] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].a[row].coeffs[pos], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[pos] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].b.coeffs[pos], decomposedMu);
        }
    }
}

void addIntegerToTrgswNtt(TrgswDft& trgswDft, Trgsw& trgsw, const Integer mu, const int pos, const YatfheParameters& param) {
    for (auto lvl = 0; lvl < trgswDft.l; lvl++) {
        auto decomposedMu = static_cast<Torus>(mu) << (param.torusBits -  (lvl + 1) * param.radixBits);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].a[row].coeffs[pos], decomposedMu);
                applyNtt(trgswDft.trlweDftSamples[lvl][row].a[row], trgsw.trlweSamples[lvl][row].a[row]);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].b.coeffs[pos], decomposedMu);
            applyNtt(trgswDft.trlweDftSamples[lvl][row].b , trgsw.trlweSamples[lvl][row].b);
        }
    }
}

void trgswAddIntegerApproxCRT(Trgsw& trgsw, const Integer mu, const YatfheParameters& param) {
    for (auto lvl = 0; lvl < trgsw.l; lvl++) {
        auto decomposedMu = multTorus(TORUS_Q, mu, param.w[lvl]);
        for (auto row = 0; row < param.k + 1; row++) {

            // add to a_lii
            if (row < param.k) {
//                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] += decomposedMu; // coeffs[0]: add mu to the constant polynomial term
                trgsw.trlweSamples[lvl][row].a[row].coeffs[0] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].a[row].coeffs[0], decomposedMu);
                continue;
            }

            // add to b_lk
//            trgsw.trlweSamples[lvl][row].b.coeffs[0] += decomposedMu;
            trgsw.trlweSamples[lvl][row].b.coeffs[0] = addTorus(TORUS_Q, trgsw.trlweSamples[lvl][row].b.coeffs[0], decomposedMu);
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
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweNoiseB, param.torusBits - param.radixBits);
}

Integer decryptTrgswNtt(const TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey) {
    const auto firstLevel = 0;
    const auto lastRow = param.k;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRoundingNtt(tmp, trgswDft.trlweDftSamples[firstLevel][lastRow], trgswKey.trlweKey);
    return roundErrorForShiftedTorus(tmp.coeffs[0], param.rlweNoiseB, param.torusBits - param.radixBits);
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
        res.coeffs[i] = roundErrorForShiftedTorus(tmp.coeffs[i], param.rlweNoiseB, param.torusBits - param.radixBits);
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
        res.coeffs[i] = roundErrorForShiftedTorus(tmp.coeffs[i], param.rlweNoiseB, param.torusBits - param.radixBits);
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

void multTrgswMPWithConst(TrgswMP& trgsw, const TrgswMP& in, const int num) {
    const auto l = trgsw.l;
    const auto k = trgsw.k;
    const auto N = trgsw.cPrime[0].b.N;
    for (auto i = 0; i < l; i++) {
        for (auto j = 0; j < k; j++) {
            for (auto z = 0; z < N; z++) {
                trgsw.cPrime[i].a[j].coeffs[z] = multTorus(TORUS_Q, in.cPrime[i].a[j].coeffs[z], num);
            }
        }
        for (auto z = 0; z < N; z++) {
            trgsw.cPrime[i].b.coeffs[z] = multTorus(TORUS_Q, in.cPrime[i].b.coeffs[z], num);
        }
        for (auto j1 = 0; j1 < k; j1++) {
            for (auto j2 = 0; j2 < k; j2++) {
                for (auto z = 0; z < N; z++) {
                    trgsw.c[i][j1].a[j2].coeffs[z] = multTorus(TORUS_Q, in.c[i][j1].a[j2].coeffs[z], num);
                }
            }
            for (auto z = 0; z < N; z++) {
                trgsw.c[i][j1].b.coeffs[z] = multTorus(TORUS_Q, in.c[i][j1].b.coeffs[z], num);
            }
        }
    }
}

void multTrgswMPWithConstNtt(TrgswMPDft& trgsw, const TrgswMPDft& in, const int num) {
    const auto l = trgsw.l;
    const auto k = trgsw.k;
    const auto N = trgsw.cPrime[0].b.N;
    const auto q = getNttHexl().GetModulus();
    // NTT of the constant polynomial [num, 0, ..., 0] evaluates to num at every
    // negacyclic root of unity, so the result is a uniform vector — no NTT needed.
    const uint64_t numMod = num >= 0 ? (uint64_t)num % q : q - (uint64_t)(-(int64_t)num) % q;
    NttPolynomial tmp{N};
    std::fill(tmp.coeffs.begin(), tmp.coeffs.end(), numMod);
    for (auto i = 0; i < l; i++) {
        for (auto j = 0; j < k; j++) {
            EltwiseMultMod(trgsw.cPrime[i].a[j].coeffs.data(), in.cPrime[i].a[j].coeffs.data(), tmp.coeffs.data(), N, q, 1);
        }
        EltwiseMultMod(trgsw.cPrime[i].b.coeffs.data(), in.cPrime[i].b.coeffs.data(), tmp.coeffs.data(), N, q, 1);
        for (auto j1 = 0; j1 < k; j1++) {
            for (auto j2 = 0; j2 < k; j2++) {
                EltwiseMultMod(trgsw.c[i][j1].a[j2].coeffs.data(), in.c[i][j1].a[j2].coeffs.data(), tmp.coeffs.data(), N, q, 1);
            }
            EltwiseMultMod(trgsw.c[i][j1].b.coeffs.data(), in.c[i][j1].b.coeffs.data(), tmp.coeffs.data(), N, q, 1);
        }
    }
}

// out = (in1 - in2) * scalar
void subMulTrgswMPNtt(TrgswMPDft& out, const TrgswMPDft& in1, const TrgswMPDft& in2, const int scalar) {
    const auto l = out.l;
    const auto k = out.k;
    const auto N = out.cPrime[0].b.N;
    const auto q = getNttHexl().GetModulus();
    const uint64_t scalarMod = scalar >= 0 ? (uint64_t)scalar % q : q - (uint64_t)(-(int64_t)scalar) % q;
    // Reuse a single scalar vector across all polynomials in the RGSW structure.
    vector<uint64_t> scalarVec(N, scalarMod);

    auto applySubMul = [&](NttPolynomial& o, const NttPolynomial& a, const NttPolynomial& b) {
        EltwiseSubMod(o.coeffs.data(), a.coeffs.data(), b.coeffs.data(), N, q);
        EltwiseMultMod(o.coeffs.data(), o.coeffs.data(), scalarVec.data(), N, q, 1);
    };
    for (auto i = 0; i < l; i++) {
        for (auto j = 0; j < k; j++) applySubMul(out.cPrime[i].a[j], in1.cPrime[i].a[j], in2.cPrime[i].a[j]);
        applySubMul(out.cPrime[i].b, in1.cPrime[i].b, in2.cPrime[i].b);
        for (auto j1 = 0; j1 < k; j1++) {
            for (auto j2 = 0; j2 < k; j2++) applySubMul(out.c[i][j1].a[j2], in1.c[i][j1].a[j2], in2.c[i][j1].a[j2]);
            applySubMul(out.c[i][j1].b, in1.c[i][j1].b, in2.c[i][j1].b);
        }
    }
}

// out = in * scalar1 + scalar2
// scalar2 is treated as a noise-free TRGSW: all a's zero, b and diagonal a's set to
// scalar2 * gadget[lvl] (uniform NTT vector) for each level.
void addMulTrgswMPWithConstNtt(TrgswMPDft& out, const TrgswMPDft& in, const int scalar1, const int scalar2) {
    const auto l = out.l;
    const auto k = out.k;
    const auto N = out.cPrime[0].b.N;
    const auto q = getNttHexl().GetModulus();

    const uint64_t s1Mod = scalar1 >= 0 ? (uint64_t)scalar1 % q : q - (uint64_t)(-(int64_t)scalar1) % q;
    NttPolynomial s1Vec{N};
    std::fill(s1Vec.coeffs.begin(), s1Vec.coeffs.end(), s1Mod);

    const uint64_t s2Mod = scalar2 >= 0 ? (uint64_t)scalar2 % q : q - (uint64_t)(-(int64_t)scalar2) % q;
    NttPolynomial addendVec{N};

    for (auto i = 0; i < l; i++) {
        const uint64_t gadgetVal = getNttGadgetRecomper(i).coeffs[0];
        const uint64_t addend = static_cast<uint64_t>((__uint128_t)s2Mod * gadgetVal % q);
        std::fill(addendVec.coeffs.begin(), addendVec.coeffs.end(), addend);

        for (auto j = 0; j < k; j++) {
            EltwiseMultMod(out.cPrime[i].a[j].coeffs.data(), in.cPrime[i].a[j].coeffs.data(), s1Vec.coeffs.data(), N, q, 1);
        }
        EltwiseMultMod(out.cPrime[i].b.coeffs.data(), in.cPrime[i].b.coeffs.data(), s1Vec.coeffs.data(), N, q, 1);
        EltwiseAddMod(out.cPrime[i].b.coeffs.data(), out.cPrime[i].b.coeffs.data(), addendVec.coeffs.data(), N, q);

        for (auto j1 = 0; j1 < k; j1++) {
            for (auto j2 = 0; j2 < k; j2++) {
                EltwiseMultMod(out.c[i][j1].a[j2].coeffs.data(), in.c[i][j1].a[j2].coeffs.data(), s1Vec.coeffs.data(), N, q, 1);
            }
            EltwiseMultMod(out.c[i][j1].b.coeffs.data(), in.c[i][j1].b.coeffs.data(), s1Vec.coeffs.data(), N, q, 1);
            EltwiseAddMod(out.c[i][j1].a[j1].coeffs.data(), out.c[i][j1].a[j1].coeffs.data(), addendVec.coeffs.data(), N, q);
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
    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(level);
    vector<TrlweDft> partials(level, TrlweDft{k, N});

    for (int lvl = 0; lvl < level; lvl++) {
        futures.emplace_back(pool.enqueue([k, N, lvl, &trgswMPInput, &decomposedTrlwe, &partials]() {
            TrlweDft inDft{k, N};
            applyNttForAB(inDft, decomposedTrlwe.trlwes[lvl]);
            const auto& c      = trgswMPInput.c[lvl];
            const auto& cPrime = trgswMPInput.cPrime[lvl];
            const auto& inA    = inDft.a;
            const auto& inB    = inDft.b;
            auto& partial      = partials[lvl];
            for (size_t i = 0; i < k; i++) {
                for (size_t i2 = 0; i2 < k; i2++) {
                    calModularInnerProductNtt(partial.a[i2], inA[i], c[i].a[i2]);
                }
                calModularInnerProductNtt(partial.b,    inA[i], c[i].b);
                calModularInnerProductNtt(partial.a[i], inB,    cPrime.a[i]);
            }
            calModularInnerProductNtt(partial.b, inB, cPrime.b);
        }));
    }
    for (auto& f : futures) f.get();

    TrlweDft tmp{k, N};
    for (int lvl = 0; lvl < level; lvl++) {
        addTrlweNtt(tmp, tmp, partials[lvl]);
    }
    applyInttForAB(output, tmp);
}

void externalProductTrgswMPNttInPlace(Trlwe& acc, const TrgswMPDft& trgswMPInput, const int level, const YatfheParameters& param) {
    const auto k = param.k;
    const auto N = param.N;
    DecomposedTrlwe decomposedTrlwe{param};
    gadgetDecomposeTrlwe(decomposedTrlwe, acc, param);

    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(level);
    vector<TrlweDft> partials(level, TrlweDft{k, N});

    for (int lvl = 0; lvl < level; lvl++) {
        futures.emplace_back(pool.enqueue([k, N, lvl, &trgswMPInput, &decomposedTrlwe, &partials]() {
            TrlweDft inDft{k, N};
            applyNttForAB(inDft, decomposedTrlwe.trlwes[lvl]);
            const auto& c      = trgswMPInput.c[lvl];
            const auto& cPrime = trgswMPInput.cPrime[lvl];
            const auto& inA    = inDft.a;
            const auto& inB    = inDft.b;
            auto& partial      = partials[lvl];
            for (size_t i = 0; i < k; i++) {
                for (size_t i2 = 0; i2 < k; i2++) {
                    calModularInnerProductNtt(partial.a[i2], inA[i], c[i].a[i2]);
                }
                calModularInnerProductNtt(partial.b,    inA[i], c[i].b);
                calModularInnerProductNtt(partial.a[i], inB,    cPrime.a[i]);
            }
            calModularInnerProductNtt(partial.b, inB, cPrime.b);
        }));
    }
    for (auto& f : futures) f.get();

    TrlweDft tmp{k, N};
    for (int lvl = 0; lvl < level; lvl++) {
        addTrlweNtt(tmp, tmp, partials[lvl]);
    }
    applyInttForAB(acc, tmp);
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
            applyNttForAB(cPrimeL, tmp);
            switchTrlweToSecretEmbeddingNtt(cL, cPrimeL, sSquare, param);
        }));
    }
    for (auto& f : futures) {
        f.wait();
    }
}

void switchTrlweToSecretEmbeddingNtt(vector<TrlweDft>& cDft, const TrlweDft& cPrimeDft, const TrlevDft& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length
    const auto N = param.N;
    auto& cPrimeADft = cPrimeDft.a;
    auto& cPrimeBDft = cPrimeDft.b;
    vector cPrimeA(K, TorusPolynomial{N});
    for (auto i = 0; i < K; i++) {
        applyIntt(cPrimeA[i], cPrimeADft[i]);
    }

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
        for (auto k2 = 0; k2 < K; k2++) {
            addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeBDft);
        }
    }
}

void switchTrlweToSecretEmbeddingNttOpt(vector<TrlweDft>& cDft, TrlweDft& cPrimeDft, const vector<vector<DecompPolynomial>>& decompA,
                                        const TrlevDft& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length
    const auto N = param.N;

    // calculate a * S^2
    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlweDfts[l];
        auto& decompL = decompA[l];
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = cDft[k1].a;
            auto& cB = cDft[k1].b;
            auto& sA = s2.a;
            auto& a = decompL[k1];
            NttPolynomial aDft{N};
            applyNtt(aDft, a);
            calModularInnerProductNtt(cPrimeDft.a[k1], aDft, getNttGadgetRecomper(l)); // calculate recomposed cPrimes'a in ntt domain
            for (auto k2 = 0; k2 < K; k2++) {
                calModularInnerProductNtt(cA[k2], aDft, sA[k2]);
            }
            calModularInnerProductNtt(cB, aDft, s2.b);
        }
    }

    // adding b
    for (auto k1 = 0; k1 < K; k1++) {
        for (auto k2 = 0; k2 < K; k2++) {
            addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeDft.b);
        }
    }
}

void switchTrlweToSecretEmbeddingNttFromDft(vector<TrlweDft>& cDft, TrlweDft& cPrimeDft,
                                            const vector<vector<NttPolynomial>>& aDft,
                                            const TrlevDft& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length

    // calculate a * S^2, using the already-NTT'd (and typically already rotated) aDft
    // supplied by the caller instead of computing NTT(decompA) here.
    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlweDfts[l];
        auto& aDftL = aDft[l];
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = cDft[k1].a;
            auto& cB = cDft[k1].b;
            auto& sA = s2.a;
            auto& a = aDftL[k1];
            calModularInnerProductNtt(cPrimeDft.a[k1], a, getNttGadgetRecomper(l)); // calculate recomposed cPrimes'a in ntt domain
            for (auto k2 = 0; k2 < K; k2++) {
                calModularInnerProductNtt(cA[k2], a, sA[k2]);
            }
            calModularInnerProductNtt(cB, a, s2.b);
        }
    }

    // adding b
    for (auto k1 = 0; k1 < K; k1++) {
        for (auto k2 = 0; k2 < K; k2++) {
            addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeDft.b);
        }
    }
}

void switchDecompTrlweToSecretEmbeddingNtt(vector<TrlweDft>& cDft, TrlweDft& cPrimeDft,
                                           const vector<vector<DecompPolynomial>>& aDecomp,
                                           const vector<DecompPolynomial>& bDecomp,
                                           const TrlevDft& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length
    const auto N = param.N;

    // calculate a * S^2, using the already-NTT'd (and typically already rotated) aDft
    // supplied by the caller instead of computing NTT(decompA) here.
    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlweDfts[l];
        auto& aL = aDecomp[l];
        auto& bL = bDecomp[l];
        NttPolynomial bDft{N};
        applyNtt(bDft, bL);
        calModularInnerProductNtt(cPrimeDft.b, bDft, getNttGadgetRecomper(l));
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = cDft[k1].a;
            auto& cB = cDft[k1].b;
            auto& sA = s2.a;
            auto& a = aL[k1];
            NttPolynomial aDft{N};
            applyNtt(aDft, a);
            calModularInnerProductNtt(cPrimeDft.a[k1], aDft, getNttGadgetRecomper(l)); // calculate recomposed cPrimes'a in ntt domain
            for (auto k2 = 0; k2 < K; k2++) {
                calModularInnerProductNtt(cA[k2], aDft, sA[k2]);
            }
            calModularInnerProductNtt(cB, aDft, s2.b);
        }
    }

    // adding b
    for (auto k1 = 0; k1 < K; k1++) {
        for (auto k2 = 0; k2 < K; k2++) {
            addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeDft.b);
        }
    }
}

void switchTrlweToSecretEmbeddingNttMix(vector<TrlweDft>& cDft, TrlweDft& cPrimeDft, const vector<vector<DecompPolynomial>>& decompA,
                                        const TorusPolynomial& cPrimeB, const TrlevDft& sSquare,
                                        const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length
    const auto N = param.N;

    // calculate a * S^2
    for (auto l = 0; l < L; l++) {
        auto& s2 = sSquare.trlweDfts[l];
        auto& decompL = decompA[l];
        for (auto k1 = 0; k1 < K; k1++) {
            auto& cA = cDft[k1].a;
            auto& cB = cDft[k1].b;
            auto& sA = s2.a;
            auto& a = decompL[k1];
            NttPolynomial aDft{N};
            applyNtt(aDft, a);
            calModularInnerProductNtt(cPrimeDft.a[k1], aDft, getNttGadgetRecomper(l)); // calculate recomposed cPrimes'a in ntt domain
            for (auto k2 = 0; k2 < K; k2++) {
                calModularInnerProductNtt(cA[k2], aDft, sA[k2]);
            }
            calModularInnerProductNtt(cB, aDft, s2.b);
        }
    }

    // adding b
    applyNtt(cPrimeDft.b, cPrimeB);
    for (auto k1 = 0; k1 < K; k1++) {
        for (auto k2 = 0; k2 < K; k2++) {
            addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeDft.b);
        }
    }
}

void switchTrlweToSecretEmbedding(vector<Trlwe>& c, const Trlwe& cPrime, const Trlev& sSquare, const YatfheParameters& param) {
    const auto K = param.k;
    const auto L = param.l; // must use full decomp length
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