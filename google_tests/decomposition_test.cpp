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
    param.l = param.torusBits / param.radixBits;
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

TEST(GADGET_DECOMP, ZERO_MEAN) {
    const auto centred = [](const int64_t d, const int bits) {
        const uint64_t Q = bits >= 64 ? 0 : (1ULL << bits);
        const uint64_t u = static_cast<uint64_t>(d) & (Q - 1);
        return u >= Q / 2 ? static_cast<int64_t>(u) - static_cast<int64_t>(Q) : static_cast<int64_t>(u);
    };
    YatfheParameters base{};
    const int tb = base.torusBits;
    const std::vector<std::pair<int, int>> gadgets{{2, tb / 2}, {4, tb / 4}, {7, 4}, {3, 5}, {3, 10}};
    const int samples = 1 << 17;
    for (const auto& [radixBits, l] : gadgets) {
        YatfheParameters param{};
        param.setRadixBits(radixBits);
        param.l = l;
        initYatfhe(param);
        const int64_t B = int64_t{1} << radixBits;
        const int dropped = tb - l * radixBits;
        const int64_t halfDelta = dropped > 0 ? int64_t{1} << (dropped - 1) : 0;
        std::vector<Torus> data(samples);
        initCoeffsViaUniformDistribution(data, TORUS_MIN, TORUS_MAX);
        std::vector<double> mean(l, 0), sq(l, 0);
        DecomposedData d{l}, dn{l};
        for (const auto x : data) {
            signedGadgetDecomposition(d, x, param);
            int64_t rec = 0;
            for (int j = 0; j < l; j++) {
                ASSERT_LE(std::llabs(static_cast<long long>(d.value[j])), B / 2);
                rec += static_cast<int64_t>(d.value[j]) * (int64_t{1} << (tb - (j + 1) * radixBits));
                mean[j] += static_cast<double>(d.value[j]);
                sq[j] += static_cast<double>(d.value[j]) * static_cast<double>(d.value[j]);
            }
            ASSERT_LE(std::llabs(static_cast<long long>(centred(static_cast<int64_t>(x) - rec, tb))), halfDelta);
            if (x != TORUS_MIN) {
                signedGadgetDecomposition(dn, static_cast<Torus>(-x), param);
                for (int j = 0; j < l; j++) ASSERT_EQ(dn.value[j], -d.value[j]);
            }
        }
        const double b = static_cast<double>(B);
        for (int j = 0; j < l; j++) {
            mean[j] /= samples;
            sq[j] /= samples;
            const int i = l - 1 - j;   // from the least significant kept level
            const double vi = (b - 2) * (b - 1) / 12.0 + b * b / (4 * (b + 1)) + (i % 2 ? -1.0 : 1.0) * std::pow(b, 1.0 - i) / (4 * (b + 1));
            // 6 standard errors of the mean; the [-B/2, B/2) digits sat at -1/2
            EXPECT_LT(std::fabs(mean[j]), 6.0 * std::sqrt(vi / samples)) << "B=2^" << radixBits << " l=" << l << " level " << j;
            if (j > 0) {   // the top digit follows the sign, not the carry chain
                EXPECT_NEAR(sq[j] / vi, 1.0, 0.03) << "B=2^" << radixBits << " l=" << l << " level " << j;
            }
        }
        // row decomposition, both unrolled layouts and the generic path, digit for digit
        const int N = 256;
        std::vector<Torus> row(data.begin(), data.begin() + N);
        std::vector<std::vector<Torus>> out(l, std::vector<Torus>(N));
        std::vector<Torus*> ptr(l);
        for (int j = 0; j < l; j++) ptr[j] = out[j].data();
        decomposeRow(ptr.data(), row.data(), N, l, param);
        for (int i = 0; i < N; i++) {
            signedGadgetDecomposition(d, row[i], param);
            for (int j = 0; j < l; j++) {
                ASSERT_EQ(out[j][i], d.value[j]);
            }
        }
    }
    // key switch digits
    YatfheParameters param{};
    param.setKsRadixBits(2);
    initYatfhe(param);
    std::vector<Torus> data(samples);
    initCoeffsViaUniformDistribution(data, LWE_MIN, LWE_MAX);
    DecomposedData d{param.ksLevel};
    std::vector<double> mean(param.ksLevel, 0);
    for (const auto x : data) {
        signedGadgetDecompositionKs(d, x, param);
        for (int j = 0; j < param.ksLevel; j++) {
            mean[j] += static_cast<double>(d.value[j]);
        }
    }
    for (int j = 0; j < param.ksLevel; j++) {
        EXPECT_LT(std::fabs(mean[j] / samples), 6.0 * std::sqrt(1.5 / samples)) << "KS level " << j;
    }
    printBanner("GADGET_DECOMP.ZERO_MEAN");
}

TEST(GADGET_DECOMP, UNSIGNED) {
    YatfheParameters param {};
    param.setRadixBits(8);
    param.l = param.torusBits / param.radixBits;
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
    param.l = param.torusBits / param.radixBits;
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
    param.l = param.torusBits / param.radixBits;
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
    param.l = param.torusBits / param.radixBits;
    param.lApprox = param.l - 1;
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

        auto z = static_cast<Torus>(longModP(recomposeTwoParts(da, bOb), TORUS_Q));
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