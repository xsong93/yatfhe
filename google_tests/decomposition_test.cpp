//
// Created by Xintong Song on 2024/3/17.
//
#include <vector>
#include <cmath>
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/crt.h"
#include "yautil/time_counter.h"
#include "yautil/initializer.h"

using namespace NttHexl;

UnsignedInteger powInt(UnsignedInteger base, UnsignedInteger exponent) {
    UnsignedInteger result = 1;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    return result;
}

TEST(GADGET_DECOMP, SIGNED) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    initYatfhe(param);
    DecomposedData decomp {param.l};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data, TORUS_MIN, TORUS_MAX);
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
        cout << "d: " << d << endl;
        printArray(decomp.value, "decomp");
        auto recons = recomposeSelf(decomp, param);
        ASSERT_EQ(d, recons);
    }
    printBanner("GADGET_DECOMP.SIGNED");
}

TEST(GADGET_DECOMP, UNSIGNED) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    initYatfhe(param);

    DecomposedData out {param.l};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data, TORUS_MIN, TORUS_MAX);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
        cout << "d: " << d << endl;
        printArray(out.value, "decomp");
        auto z = recomposeSelf(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("GADGET_DECOMP.UNSIGNED");
}

TEST(GADGET_DECOMP, DECOMPOSE_OVER_B_SINGLE_STAGE) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    initYatfhe(param);

    Torus mult = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    vector<Torus> rhs(param.l);
    decomposeOverB(rhs, mult, param.torusBits, param.radixBits);
    printArray(rhs, to_string(mult) + " decomposeOverB");
    DecomposedData decomp {param.l};
    int ti = 0;
    while (ti++ < 10) {
        Torus data = genIntUniformDist(TORUS_MIN, TORUS_MAX);
        signedGadgetDecomposition(decomp, data, param); // both correct
        // gadgetDecompose(decomp, data, param); // both correct
        cout << endl << "iter: " << ti << ", " << "in: " << data << endl;
        printArray(decomp.value, "decomp");
        auto out = recomposeTwoParts(decomp, rhs);
        cout << "out: " << out << ", " << "data * mult: " << data * mult << endl;
        ASSERT_EQ(out, data * mult);
    }
    printBanner("GADGET_DECOMP.DECOMPOSE_OVER_B_SINGLE_STAGE");
}

TEST(GADGET_DECOMP, DECOMPOSE_OVER_B_MULTI_STAGE) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    initYatfhe(param);
    auto level = param.l;

    Torus mult = 5;

    // first decomp
    std::vector<Torus> rhs(level);
    decomposeOverB(rhs, mult, param.torusBits, param.radixBits);
    printArray(rhs, to_string(mult) + " decomposeOverB");

    // second decomp
    vector<DecomposedData> decompRhs {rhs.size(), DecomposedData(level)};
    for (auto i = 0; i < rhs.size(); i++) {
        gadgetDecompose(decompRhs[i], rhs[i], param);
        printArray(decompRhs[i].value, "rhs level " + to_string(i));
    }

    // recomp second decomposed data test
    vector<Integer> recompRhs(rhs.size());
    for (auto i = 0; i < decompRhs.size(); i++) {
        recompRhs[i] = recomposeSelf(decompRhs[i], param);
        ASSERT_EQ(recompRhs[i], rhs[i]);
    }
    printArray(recompRhs, "recompRhs");

    DecomposedData decompL1 {level};
    DecomposedData recompL1 {level};
//    Torus data = genIntUniformDist(TorusMin, TorusMax);
    Torus data = TORUS_MAX - 3;

//  signedGadgetDecomposition(decomp, data, param); // both correct
    gadgetDecompose(decompL1, data, param); // both correct
    cout << "in: " << data << endl;
    printArray(decompL1.value, "decomp");

    // recomp first level
    recomposeFirstHalf(recompL1, decompL1, decompRhs);
    printArray(recompL1.value, "recompL1");

    // recomp second level
    auto out = recomposeSelf(recompL1, param); // equivalent to recomposeTwoParts(recompL1, decompOneOverR)
    cout << "out: " << out << ", " << "data * mult: " << data * mult << endl;

    for (auto i = 0; i < level; i++) {
        cout << decompL1.value[i] * rhs[i] << ", "
             << ((decompL1.value[i] * rhs[i]) >> (param.torusBits - level * (i+1))) << endl;
    }

    ASSERT_EQ(out, data * mult);

    printBanner("GADGET_DECOMP.DECOMPOSE_OVER_B_MULTI_STAGE");
}

TEST(GADGET_DECOMP, DECOMPOSE_TRLWE) {
    YatfheParameters param {};
    param.torusBase = 8;
    param.setRadixBits(8);
    param.l = 4;
    initYatfhe(param);
    auto level = param.l;

    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};
    TrlweDft recompDft {param.k, param.N};

    int plain = 1;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TrlweKey trlweKey {param};
    genTrlweKey(trlweKey);
//    Trlwe in2 {param.k, param.N};
//    symEncTrlweSingleSample(in2, trlweKey, mu);
    int pos = 0;
    symEncTrlweSingleSampleNtt(in, inDft, trlweKey, mu, 0);

    // decompose
    DecomposedTrlwe out {param, level};
    gadgetDecomposeTrlwe(out, in, param);

    // recompose original
    recomposeTrlwe(recomp, out, param);

    IntPolynomial dec1 {param.N};
    symDecTrlweToInt(dec1, in, trlweKey, param.torusBase);
    for (auto j = 0; j < param.N; j++) {
        if (j == pos) {
            ASSERT_EQ(plain, dec1.coeffs[pos]);
            continue;
        }
        ASSERT_EQ(0, dec1.coeffs[j]);
    }

    printBanner("GADGET_DECOMP.DECOMPOSE_TRLWE");
}

TEST(GADGET_DECOMP, ADD_SUB) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    auto level = param.l;

    DecomposedData da {level};
    DecomposedData db {level};
    DecomposedData dr {level};
    int a;
    int b;
    int t = 5000;
    while (t-- > 0) {
        a = genIntUniformDist(INT32_MIN, INT32_MAX);
        b = genIntUniformDist(INT32_MIN, INT32_MAX);
        signedGadgetDecomposition(da, a, param);
        signedGadgetDecomposition(db, b, param);
        for (auto i = 0; i < level; i++) {
            dr.value[i] = da.value[i] * da.sign + db.value[i] * db.sign;
        }
        auto z = recomposeSelf(dr, param);
        cout << "a + b: a: " << a     << ", "
             << "b: "        << b     << ", "
             << "decomp: "   << z     << ", "
             << "ori: "      << a + b << endl;
        for (auto i = 0; i < level; i++) {
            dr.value[i] = da.value[i] * da.sign - db.value[i] * db.sign;
        }
        z = recomposeSelf(dr, param);
        cout << "a - b: a: " << a     << ", "
             << "b: "        << b     << ", "
             << "decomp: "   << z     << ", "
             << "ori: "      << a - b << endl;
    }
    printBanner("GADGET_DECOMP.ADD_SUB");
}


// Actually, identical logic with DecomposeOverB test.
TEST(GADGET_DECOMP, MULT) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = 4;
    param.lApprox = 3;
    initYatfhe(param);
    auto level = param.l;

    DecomposedData da {level};
    DecomposedData db {level};
    DecomposedData dr {level};
    Torus a;
    Torus b;
    int t = 2000;
    while (t-- > 0) {
        std::cout << "iter: " << t << endl;
        a = genIntUniformDist(TORUS_MIN, TORUS_MAX);
        b = genIntUniformDist(TORUS_MIN, TORUS_MAX);

        signedGadgetDecomposition(da, a, param);
        printArray(da.value, "da");

        std::vector<Torus> bOb(level);
        decomposeOverB(bOb, b, param.torusBits, param.radixBits);
        printArray(bOb, "bOb");

        auto z = recomposeTwoParts(da, bOb);
        DecomposedData dz {level};
        signedGadgetDecomposition(dz, z, param);
        printArray(dz.value, "dz");

        auto r = multTorus(TORUS_Q, a, b);
        DecomposedData dt {level};
        signedGadgetDecomposition(dt, r, param);
        printArray(dt.value, "dt");

        cout << "a * b: a: " << a     << ", "
             << "b: "        << b     << ", "
             << "decomp: "   << z     << ", "
             << "ori: "      << r     << endl;
        ASSERT_EQ(z, r);
    }
    printBanner("GADGET_DECOMP.MULT");
}