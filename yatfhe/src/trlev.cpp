//
// Created by Xintong Song on 2024/4/22.
//
#include <vector>
#include "yatfhe/trlev.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/numeric.h"

using namespace std;

void encTrlevSingleSample(Trlev& output, const TrlweKey& trlweKey, const Torus input, const int pos, const YatfheParameters& param) {
    const auto l = output.l;
    for (auto i = 0; i < l; i++) {
        auto inOverR = input << (param.torusBits - (i + 1) * param.radixBits);
        symEncTrlweSingleSample(output.trlwes[i], trlweKey, inOverR, pos);
    }
}

void encTrlevSingleSampleMonomial(Trlev& output, const TrlweKey& trlweKey, const Torus input, const int monomialIndex,
                                   const YatfheParameters& param) {
    const auto l = output.l;
    const auto N = param.N;
    for (auto lvl = 0; lvl < l; lvl++) {
        const Torus inOverR = input << (param.torusBits - (lvl + 1) * param.radixBits);
        Trlwe& ct = output.trlwes[lvl];

        // Initialize b as the RLWE plaintext polynomial v = inOverR * X^{monomialIndex} (+ noise),
        // and initialize a uniformly. Then apply symEncTrlwe to add the a*s term.
        for (auto j = 0; j < N; j++) {
            const Torus msg = (j == monomialIndex) ? inOverR : 0;
            ct.b.coeffs[j] = addTUniformNoise(msg, trlweKey.errorB, TORUS_Q);
        }
        for (auto k = 0; k < ct.k; k++) {
            initCoeffsViaUniformDistribution(ct.a[k].coeffs, TORUS_MIN, TORUS_MAX);
        }
        for (auto k = 0; k < ct.k; k++) {
            multTorusPolynomialAcc(ct.b, ct.a[k], trlweKey.s[k]);
        }
    }
}

void encTrlevMultiSample(Trlev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param) {
    const auto l = output.l;
    vector<Torus> inputsOverR(inputs.N);
    for (auto i = 0; i < l; i++) {
        for (auto j = 0; j < inputs.N; j++) {
            inputsOverR[j] = inputs.coeffs[j] << (param.torusBits - (i + 1) * param.radixBits);
        }
        symEncTrlweMultiSample(output.trlwes[i], trlweKey, inputsOverR);
    }
}

// encrypt +-s*m, sign is determined by @param{isPos}
// (a-x, as+e)
void symEncTrlevWithKey(Trlev& output, const TrlweKey& trlweKey, const vector<TorusPolynomial>& inputs, const bool isPos,
                        const YatfheParameters& param) {
    encTrlevSingleSample(output, trlweKey, 0, 0, param);
    for (auto i = 0; i < param.l; i++) {
        for (auto j = 0; j < param.k; j++) {
            TorusPolynomial sXm{param.N};
            for (auto z = 0; z < param.N; z++) {
                sXm.coeffs[z] = inputs[j].coeffs[z] << (param.torusBits - (i + 1) * param.radixBits);
            }
            if (isPos) {
                subTorusPolynomial(output.trlwes[i].a[j], output.trlwes[i].a[j], sXm);
            } else {
                addTorusPolynomial(output.trlwes[i].a[j], output.trlwes[i].a[j], sXm);
            }
        }
    }
}

// encrypt +-s*m to ntt domain, sign is determined by @param{isPos}
// (a-x, as+e)
void symEncTrlevWithKeyNtt(TrlevDft& output, const TrlweKey& trlweKey, const vector<TorusPolynomial>& inputs,
                           const bool isPos, const YatfheParameters& param) {
    Trlev s2(param);
    encTrlevSingleSample(s2, trlweKey, 0, 0, param);
    for (auto i = 0; i < param.l; i++) {
        NttHexl::applyNtt(output.trlweDfts[i].b, s2.trlwes[i].b);
        for (auto j = 0; j < param.k; j++) {
            TorusPolynomial sXm{param.N};
            for (auto z = 0; z < param.N; z++) {
                sXm.coeffs[z] = inputs[j].coeffs[z] << (param.torusBits - (i + 1) * param.radixBits);
            }
            if (isPos) {
                subTorusPolynomial(s2.trlwes[i].a[j], s2.trlwes[i].a[j], sXm);
            } else {
                addTorusPolynomial(s2.trlwes[i].a[j], s2.trlwes[i].a[j], sXm);
            }
            NttHexl::applyNtt(output.trlweDfts[i].a[j], s2.trlwes[i].a[j]);
        }
    }
}

void decTrlev(TorusPolynomial& output, const Trlev& input, const TrlweKey& trlweKey, const YatfheParameters& param) {
    const auto firstLevel = 0;
    TorusPolynomial tmp {param.N};
    symDecTrlweWoRounding(tmp, input.trlwes[firstLevel], trlweKey);
    for (auto i = 0; i < param.N; i++) {
        output.coeffs[i] = roundErrorForShiftedTorus(tmp.coeffs[i], param.rlweNoiseB, param.torusBits - param.radixBits);
    }
}

void multTrlevWithConst(Trlwe& output, const Trlev& input, const Torus num, const YatfheParameters& param) {
    auto N = output.b.N;
    auto k = output.k;
    DecomposedData d {input.l};
    signedGadgetDecomposition(d, num, param);
    for (auto j = 0; j < N; j++) {
        for (auto r = 0; r < k + 1; r++) {
            auto& curr = (r < k) ? output.a[r] : output.b;
            for (auto l1 = 0; l1 < d.l; l1++) {
                auto& currTglev = (r < k) ? input.trlwes[l1].a[r] : input.trlwes[l1].b;
//                curr.coeffs[j] += currTglev.coeffs[j] * d.value[l1] * d.sign;
                curr.coeffs[j] = addTorus(TORUS_Q, curr.coeffs[j], multTorus(TORUS_Q, currTglev.coeffs[j], d.value[l1] * d.sign));
            }
        }
    }
}

void multDecomposedTrlevWithConst(Trlwe& output, const Trlev& input, const Torus num, const YatfheParameters& param) {
    const auto N = output.b.N;
    const auto k = output.k;
    const auto lvl0 = input.l;
    const auto lvl1 = param.l2;
    DecomposedData lhs {lvl0};
    signedGadgetDecomposition(lhs, num, param);
    std::vector<DecomposedTrlwe> decompTglev(lvl1, DecomposedTrlwe(param));
    for (auto i = 0; i < lvl1; i++) {
        gadgetDecomposeTrlwe(decompTglev[i], input.trlwes[i], param);
    }

    // first recomp
    DecomposedTrlwe recomp1 {param};
    for (auto l0 = 0; l0 < lvl0; l0++) {
        for (auto x1 = 0; x1 < N; x1++) {
            for (auto x2 = 0; x2 < k + 1; x2++) {
                auto& curr = (x2 < k) ? recomp1.trlwes[l0].a[x2] : recomp1.trlwes[l0].b;
                for (auto l1 = 0; l1 < lvl1; l1++) {
                    auto& currDecomp = (x2 < k) ? decompTglev[l1].trlwes[l0].a[x2] : decompTglev[l1].trlwes[l0].b;
                    curr.coeffs[x1] += (lhs.value[l1] * lhs.sign) * (currDecomp.coeffs[x1]);
                }
            }
        }
    }

    recomposeTrlwe(output, recomp1, param);
}

void rotateTrlev(Trlev& trglev, const int rot, const YatfheParameters& param) {
    if (rot % (param.N * 2) == 0) {
        return;
    }
    Trlwe tmp{param.k, param.N};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        tmp = trglev.trlwes[lvl];
        rotateTrlwe(trglev.trlwes[lvl], tmp, rot);
    }
}