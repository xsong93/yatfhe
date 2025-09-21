//
// Created by Xintong Song on 2024/4/22.
//
#include <vector>
#include "yatfhe/trlev.h"

using namespace std;

void encTrlevSingleSample(Trlev& output, const TrlweKey& trlweKey, const Torus input, const YatfheParameters& param) {
    const auto l = output.l;
    for (auto i = 0; i < l; i++) {
        auto inOverR = input << (param.torusBits - (i + 1) * param.radixBits);
        symEncTrlweSingleSample(output.trlwes[i], trlweKey, inOverR);
    }
}

void encTrlevMultiSample(Trlev& output, const TrlweKey& trlweKey, const TorusPolynomial& inputs, const YatfheParameters& param) {
    const auto l = output.l;
    vector<Torus> inputsOverR(inputs.N);
    for (auto i = 0; i < l; i++) {
        for (auto j = 0; j < inputs.N; j++) {
            inputsOverR[j] = inputs.coeffs[j] << (param.torusBits - (i + 1) * param.radixBits);;
        }
        symEncTrlweMultiSample(output.trlwes[i], trlweKey, inputsOverR);
    }
}

void multTrlevWithConst(Trlwe& output, const Trlev& input, const Integer num, const YatfheParameters& param) {
    auto N = output.b.N;
    auto k = output.k;
    DecomposedData d {input.l};
    gadgetDecompose(d, num, param);
    for (auto j = 0; j < N; j++) {
        for (auto r = 0; r < k + 1; r++) {
            auto& curr = (r < k) ? output.a[r] : output.b;
            for (auto l1 = 0; l1 < d.l; l1++) {
                auto& currTglev = (r < k) ? input.trlwes[l1].a[r] : input.trlwes[l1].b;
//                curr.coeffs[j] += currTglev.coeffs[j] * d.value[l1] * d.sign;
                curr.coeffs[j] = addTorus(curr.coeffs[j], multTorus(currTglev.coeffs[j], d.value[l1] * d.sign));
            }
        }
    }
}

void multDecomposedTrlevWithConst(Trlwe& output, const Trlev& input, const Integer num, const YatfheParameters& param) {
    const auto N = output.b.N;
    const auto k = output.k;
    const auto lvl0 = input.l;
    const auto lvl1 = param.l2;
    DecomposedData lhs {lvl0};
    gadgetDecompose(lhs, num, param);
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