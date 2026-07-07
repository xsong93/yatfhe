//
// Created by Xintong Song on 2023/12/8.
//
#include <iostream>
#include "yautil/control_helper.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/numeric.h"

using namespace std;

void genTlweKey(TlweKey& key) {
    for (auto i = 0; i < key.n; i++) {
#ifdef TERNRY
        key.s[i] = ternaryDistrib(rng);
#else
        key.s[i] = binaryDistrib(rng);
#endif
    }
}

// b = aj * sj + u + e
void symEncTlwe(Tlwe& tlweSample, const Torus message, const TlweKey& key) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        tlweSample.a[i] = uniformTorusDistrib(LWE_MIN, LWE_MAX)(rng);
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(tlweSample.a[i]) * key.s[i];
        }
    }
    Torus muE = addTUniformNoise(message, key.errorB, LWE_Q);
    tlweSample.b = addTorus(LWE_Q, muE, static_cast<Torus>(longModP(tmp, LWE_Q)));
}

Torus symDecTlweToTorus(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, LWE_Q));
    return roundTorusGeneralError(subTorus(LWE_Q, in.b, aXs), torusBase, LWE_Q);
}

Integer symDecTlweToInt(Tlwe& in, const TlweKey& key, const int torusBase) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, LWE_Q));
    return modSwitchFromTorusGeneral(roundTorusGeneralError(subTorus(LWE_Q, in.b, aXs), torusBase, LWE_Q), torusBase, LWE_Q);
}

Torus calTlweError(Tlwe& in, const TlweKey& key, const Integer mu) {
    int64_t tmp = 0;
    for (auto i = 0; i < key.n; i++) {
        if (key.s[i] != 0) {
            tmp += static_cast<int64_t>(in.a[i]) * key.s[i];
        }
    }
    auto aXs = static_cast<Torus>(longModP(tmp, LWE_Q));
    return subTorus(LWE_Q, in.b, aXs) - modSwitchToTorusGeneral(mu, MESSAGE_P, LWE_Q);
}

void rescaleTlweFromTorus(Tlwe& output, const Tlwe& input) {
    output.b = static_cast<int32_t>(modSwitchFromTorusGeneral(input.b, LWE_Q, TORUS_Q));
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = static_cast<int32_t>(modSwitchFromTorusGeneral(input.a[i], LWE_Q, TORUS_Q));
    }
}

void rescaleTlweToNewMod(ScaledTlwe& output, const Tlwe& input) {
    const auto newMod = output.mod;
    output.b = static_cast<int32_t>(modSwitchFromTorusGeneral(input.b, newMod, LWE_Q));
    for (auto i = 0; i < input.n; i++) {
        output.a[i] = static_cast<int32_t>(modSwitchFromTorusGeneral(input.a[i], newMod, LWE_Q));
    }
}

void addTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = addTorus(LWE_Q, input1.a[i], input2.a[i]);
    }
    output.b = addTorus(LWE_Q, input1.b, input2.b);
}

void subTlwe(Tlwe& output, const Tlwe& input1, const Tlwe& input2) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = subTorus(LWE_Q, input1.a[i], input2.a[i]);
    }
    output.b = subTorus(LWE_Q, input1.b, input2.b);
}

// output -= input
void subTlweInPlace(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = subTorus(LWE_Q, output.a[i], input.a[i]);
    }
    output.b = subTorus(LWE_Q, output.b, input.b);
}

void multTlwe(Tlwe& output, const int scalar) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = multTorus(LWE_Q, output.a[i], scalar);
    }
    output.b = multTorus(LWE_Q, output.b, scalar);
}

void copyTlwe(Tlwe& output, const Tlwe& input) {
    for (auto i = 0; i < output.n; i++) {
        output.a[i] = input.a[i];
    }
    output.b = input.b;
}

void resetTlweToZero(Tlwe& tlwe) {
    tlwe.b = 0;
    for (auto& ai : tlwe.a) {
        ai = 0;
    }
}

void inverseTlwe(Tlwe& tlwe) {
    tlwe.b = modMulQ(tlwe.b, static_cast<Torus>(-1), LWE_Q);
    for (auto i = 0; i < tlwe.n; i++) {
        tlwe.a[i] = modMulQ(tlwe.a[i], static_cast<Torus>(-1), LWE_Q);
    }
}
