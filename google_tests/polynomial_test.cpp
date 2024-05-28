//
// Created by Xintong Song on 2024/2/21.
//

#include "gtest/gtest.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/tlwe.h"
#include "yautil/tool.h"

TEST(PolynomialTest, PolynomialRounding) {
    YatfheParameters param {};
//    param.torusBase = 8;
//    param.N = 1024;
//    param.n = 4;
    TorusPolynomial v {param.N};
    std::vector<Integer> d(param.N);
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    printArray(v.coeffs, "v");
    for (int i = 0; i < v.N; i++) {
        d[i] = modSwitchFromTorus32(v.coeffs[i], param.torusBase);
    }
    printArray(d, "p");

    Integer plain = 3;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TlweKey tlweKey {param.n, 0};
    lweKeyGen(tlweKey);

    Tlwe ct {param.n};
    symEncTlweSample(ct, mu, tlweKey);

    ScaledTlwe scaledCt {2 * param.N, param.n};
    rescaleTlweFromTorus32(scaledCt, ct);
    printTlweAB(ct, "ct");
    printTlweAB(scaledCt, "scaledCt");

    TorusPolynomial rpT {param.N};
    int rot = scaledCt.b;
    for (auto i = 0; i < scaledCt.n; i++) {
        rot = (rot - scaledCt.a[i] * tlweKey.s[i]) % (2 * param.N);
    }
    cout << "u*: " << rot << endl;
    torusPolynomialRotate(rpT, -rot, v);
    IntPolynomial res {param.N};
    torusPolyToIntPoly(res, rpT, param.torusBase);
    cout << "p: " << intModP(plain, param.torusBase) << endl;
    printArray(res.coeffs, "res");
    printBanner("PolynomialRounding");
}