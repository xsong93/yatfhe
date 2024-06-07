//
// Created by Xintong Song on 2024/3/17.
//
#include <vector>
#include <cmath>
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/gadget_decomposition.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/ntt.h"
#include "yautil/time_counter.h"
#include "yautil/initializer.h"

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

TEST(DecompositionTest, SignedDecompTest) {
    YatfheParameters param {};
    DecomposedData decomp {param.ksLevel};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data);
//    data = {-10,-9,9,10};
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
        printf("in: %d, ", d);
        printArray(decomp.value, "decomp");
        auto recons = selfRecompose(decomp, param);
        ASSERT_EQ(d, recons);
    }
    printBanner("SignedDecomp");
}

TEST(DecompositionTest, DecomposeTest) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;
    DecomposedData out {param.ksLevel};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
        printf("in: %d, ", d);
        printArray(out.value, "decomp");
        auto z = selfRecompose(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("Decompose");
}

TEST(DecompositionTest, DecomposeOverBSingleStage) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;
    Torus mult = genIntUniformDist(TorusMin, TorusMax);
    std::vector<Integer> rhs(param.ksLevel);
    decomposeOverB(rhs, mult, param);
    printArray(rhs, to_string(mult) + " decomposeOverB");
    DecomposedData decomp {param.ksLevel};
    int ti = 0;
    while (ti++ < 10) {
        Torus data = genIntUniformDist(TorusMin, TorusMax);
//        signedGadgetDecomposition(decomp, data, param); // both correct
        gadgetDecompose(decomp, data, param); // both correct
        printf("iter: %d, in: %d, ", ti, data);
        printArray(decomp.value, "decomp");
        auto out = recomposeTwoParts(decomp, rhs);
        printf("out: %d, data * mult: %d\n\n", out, data * mult);
        ASSERT_EQ(out, data * mult);
    }
    printBanner("DecomposeOverBSingleStage");
}

TEST(DecompositionTest, DecomposeOverBMultiStages) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;

    Torus mult = 5;

    // first decomp
    std::vector<Integer> rhs(param.ksLevel);
    decomposeOverB(rhs, mult, param);
    printArray(rhs, to_string(mult) + " decomposeOverB");

    // second decomp
    vector<DecomposedData> decompRhs {rhs.size(), DecomposedData(param.ksLevel)};
    for (auto i = 0; i < rhs.size(); i++) {
        gadgetDecompose(decompRhs[i], rhs[i], param);
        printArray(decompRhs[i].value, "rhs level " + to_string(i));
    }

    // recomp second decomposed data test
    vector<Integer> recompRhs(rhs.size());
    for (auto i = 0; i < decompRhs.size(); i++) {
        recompRhs[i] = selfRecompose(decompRhs[i], param);
        ASSERT_EQ(recompRhs[i], rhs[i]);
    }
    printArray(recompRhs, "recompRhs");

    DecomposedData decompL1 {param.ksLevel};
    DecomposedData recompL1 {param.ksLevel};
//    Torus data = genIntUniformDist(TorusMin, TorusMax);
    Torus data = TorusMax - 3;

//  signedGadgetDecomposition(decomp, data, param); // both correct
    gadgetDecompose(decompL1, data, param); // both correct
    printf("in: %d, ", data);
    printArray(decompL1.value, "decomp");

    // recomp first level
    recomposeFirstHalf(recompL1, decompL1, decompRhs);
    printArray(recompL1.value, "recompL1");

    // recomp second level
    auto out = selfRecompose(recompL1, param); // equivalent to recomposeTwoParts(recompL1, decompOneOverR)
    printf("out = %d, data * mult = %d\n", out, data * mult);

    for (auto i = 0; i < param.ksLevel; i++) {
        cout << decompL1.value[i] * rhs[i] << ", "<< ((decompL1.value[i] * rhs[i]) >> (32 - 4* (i+1))) << endl;
    }

    ASSERT_EQ(out, data * mult);

    printBanner("DecomposeOverBMultiStages");
}

TEST(DecompositionTest, DecomposeTrlweTest) {
    YatfheParameters param {};
    yatfheInit(param);

    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};
    TrlweDft recompDft {param.k, param.N};

    Integer plain = 1;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    trlweKeyGen(trlweKey);
    symEncTrlweSingleSampleNtt(in, inDft, trlweKey, mu);
    Trlwe intt {param.k, param.N};
    applyInttForAB(intt, inDft);

    // decompose
    DecomposedTrlwe out {param};
    DecomposedTrlweDft outDft {param, param.lDft}; // lDft > l, optimize?
    gadgetDecomposeTrlwe(out, in, param);
    gadgetDecomposeTrlweNtt(outDft, inDft, param);
    printDecomposedTrlweAB(out, "out");
    printDecomposedTrlweNttAB(outDft, "outNtt");

    // recompose original
    recomposeTrlwe(recomp, out, param);
    printTrlweAB(in, "original");
    printTrlweAB(recomp, "recomp");
    IntPolynomial dec1 {param.N};
    symDecTrlweToInt(dec1, recomp, trlweKey, param.torusBase);
    for (auto j = 0; j < in.b.N; j++) {
        ASSERT_EQ(dec1.coeffs[j], plain);
    }

    // recompose ntt
    recomposeTrlweNtt(recompDft, outDft, param);
    printTrlweDftAB(inDft, "inDft");
    printTrlweDftAB(recompDft, "recompDft");
    symDecTrlweToIntNtt(dec1, recompDft, trlweKey, param.torusBase);
    for (auto j = 0; j < in.b.N; j++) {
        ASSERT_EQ(dec1.coeffs[j], plain);
    }
    printBanner("DecomposeTrlweTest");
}

TEST(DecompositionTest, DecomposedAddSub) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;
    DecomposedData da {param.ksLevel};
    DecomposedData db {param.ksLevel};
    DecomposedData dr {param.ksLevel};
    int a;
    int b;
    int t = 5000;
    while (t-- > 0) {
        a = genIntUniformDist(INT32_MIN, INT32_MAX);
        b = genIntUniformDist(INT32_MIN, INT32_MAX);
        gadgetDecompose(da, a, param);
        gadgetDecompose(db, b, param);
        for (auto i = 0; i < param.ksLevel; i++) {
            dr.value[i] = da.value[i] * da.sign + db.value[i] * db.sign;
        }
        auto z = selfRecompose(dr, param);
        printf("a + b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a + b);
        for (auto i = 0; i < param.ksLevel; i++) {
            dr.value[i] = da.value[i] * da.sign - db.value[i] * db.sign;
        }
        z = selfRecompose(dr, param);
        printf("a - b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a - b);
    }
    printBanner("DecomposedAddSub");
}


// Actually, identical logic with DecomposeOverB test.
TEST(DecompositionTest, DecomposedMult) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;
    DecomposedData da {param.ksLevel};
    DecomposedData db {param.ksLevel};
    DecomposedData dr {param.ksLevel};
    Torus a;
    Torus b;
    int t = 2000;
    while (t-- > 0) {
        std::cout << "iter: " << t << endl;
        a = genIntUniformDist(TorusMin, TorusMax);
        b = genIntUniformDist(TorusMin, TorusMax);

        gadgetDecompose(da, a, param);
        printArray(da.value, "da");

        std::vector<Integer> bOb(param.ksLevel);
        decomposeOverB(bOb, b, param);
        printArray(bOb, "bOb");

        auto z = recomposeTwoParts(da, bOb);
        DecomposedData dz {param.ksLevel};
        gadgetDecompose(dz, z, param);
        printArray(dz.value, "dz");

//        auto r = (int)(((int64_t)a * (int64_t)b) % (param.q));
        auto r = a * b;
        DecomposedData dt {param.ksLevel};
        gadgetDecompose(dt, r, param);
        printArray(dt.value, "dt");

        printf("a * b: a: %d, b: %d. decomp: %d, ori: %d\n\n\n", a, b, z, r);
        ASSERT_EQ(z, r);
    }
    printBanner("DecomposedMult");
}

TEST(DecompositionTest, DecompNttOrderTest) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.l = 8;
//    param.k = 2;
    yatfheInit(param);

    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    trlweKeyGen(trlweKey);
    symEncTrlweSingleSampleNtt(in, inDft, trlweKey, mu);

//    TorusPolynomial helper {param.N, 1};
//    LagrangePolynomial helperDft {param.N};
//    applyNtt(helperDft, helper);
//    printArray(helperDft.coeffs, "helperDft");
//    vector<DecomposedDataDft> decomposedHelperDft;
//    for (auto i = 0; i < helperDft.N; i++) {
//        DecomposedDataDft ddft {param.lDft};
//        gadgetDecomposeNtt(ddft, helperDft.coeffs[i], param);
//        decomposedHelperDft.push_back(ddft);
//    }


    // NTT -> decomp -> recomp -> INTT
    DecomposedTrlwe out1 {param};
    DecomposedTrlweDft outDft1 {param, param.lDft};
    TrlweDft dft1 {param.k, param.N};
    TrlweDft recompDft1 {param.k, param.N};
    Trlwe intt1 {param.k, param.N};
    COUNT_TIME("NTT", applyNttForAB(dft1, in);) // NTT
    COUNT_TIME("decomp", gadgetDecomposeTrlweNtt(outDft1, dft1, param);) //decomp
    COUNT_TIME("recomp", recomposeTrlweNtt(recompDft1, outDft1, param);) // recompose
    COUNT_TIME("INTT", applyInttForAB(intt1, recompDft1);) // intt

    TorusPolynomial dec1 {param.N};
    symDecTrlweToInt(dec1, intt1, trlweKey, param.torusBase);
    printArray(dec1.coeffs, "dec1");
    for (auto j = 0; j < in.b.N; j++) {
        for (auto i = 0; i < in.k; i++) {
            ASSERT_EQ(in.a[i].coeffs[j], intt1.a[i].coeffs[j]);
        }
        ASSERT_EQ(in.b.coeffs[j], intt1.b.coeffs[j]);
    }
    printBanner("NTT -> decomp -> recomp -> INTT");


//    // decomp -> NTT -> recomp -> INTT
//    // todo
//    DecomposedTrlwe out2 {param};
//    DecomposedTrlweDft outDft2 {param, param.l};
//    TrlweDft recompDft2 {param.k, param.N};
//    Trlwe intt2 {param.k, param.N};
//
//    std::vector<Integer> oneOverB(param.l);
//    decomposeOverB(oneOverB, 1, param);
//    printArray(oneOverB, "1");
//    vector<TorusPolynomial> helper (param.l, TorusPolynomial(param.N, 0));
//    vector<LagrangePolynomial> helperDft;
//    for (auto i = 0; i < helper.size(); i++) {
//        for (auto j = 0; j < param.N; j++) {
//            helper[i].coeffs[j] = oneOverB[i];
//        }
//        LagrangePolynomial tmp (param.N);
//        applyNtt(tmp, helper[i]);
//        helperDft.push_back(tmp);
//    }
//    for (auto i = 0; i < helperDft.size(); i++) {
//        printArray(helperDft[i].coeffs, "helperDft l" + to_string(i));
//    }
//
//    COUNT_TIME("decomp", gadgetDecomposeTrlwe(out2, in, param);) // decompose
//    COUNT_TIME("NTT",
//               for (auto i = 0; i < out2.l; i++) {
//                   applyNttForAB(outDft2.rlweDfts[i], out2.rlwes[i]); // NTT
//               })
//    COUNT_TIME("recomp", {
//        auto l = out2.l;
//        auto k = param.k;
//        for (auto lvl = 0; lvl < l; lvl++) {
//            for (auto row = 0; row < k + 1; row++) {
//                auto& currIn = (row < k) ? outDft2.rlweDfts[lvl].a[row] : outDft2.rlweDfts[lvl].b;
//                auto& currOut = (row < k) ? recompDft2.a[row] : recompDft2.b;
//                for (auto j = 0; j < param.N; j++) {
//                    currOut.coeffs[j] = modAdd(currOut.coeffs[j], modMul(currIn.coeffs[j], helperDft[lvl].coeffs[j]));
//                }
//            }
//        }
//    }) // recompose
//    COUNT_TIME("intt", applyInttForAB(intt2, recompDft2);) // intt
//    TorusPolynomial dec2 {param.N};
//    symDecTrlweToInt(dec2, intt2, trlweKey, param.torusBase);
//    printArray(dec2.coeffs, "dec");
//    for (auto j = 0; j < in.b.N; j++) {
//        for (auto i = 0 ; i < in.k; i++) {
//            ASSERT_EQ(in.a[i].coeffs[j], intt2.a[i].coeffs[j]);
//        }
//        ASSERT_EQ(in.b.coeffs[j], intt2.b.coeffs[j]);
//    }
//    printBanner("decomp -> NTT -> recomp -> INTT");
}