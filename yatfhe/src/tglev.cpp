//
// Created by Xintong Song on 2024/6/12.
//
#include "yatfhe/tglev.h"
#include "yatfhe/gadget_decomposition.h"

using namespace std;

void tglevEnc(Tglev& output, const TlweKey& tlweKey, Torus input, const YatfheParameters& param) {
    const auto l = output.l;
    for (auto i = 0; i < l; i++) {
        auto inOverR = input << (param.torusBits - (i + 1) * param.radixBits);
        symEncTlweSample(output.tlwes[i], inOverR, tlweKey);
    }
}

void tglevMultConst(Tlwe& output, const Tglev& input, Integer num, const YatfheParameters& param) {
    DecomposedData d {input.l};
    gadgetDecompose(d, num, param);
    for (auto j = 0; j <= output.n; j++) {
        auto& curr = (j < output.n) ? output.a[j] : output.b;
        for (auto l1 = 0; l1 < d.l; l1++) {
            auto& currTglev = (j < output.n) ? input.tlwes[l1].a[j] : input.tlwes[l1].b;
            curr += currTglev * d.value[l1] * d.sign;
        }
    }
}