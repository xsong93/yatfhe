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

TEST(DecompositionTest, SignedDecompTest) {
    YatfheParameters param {};
    DecomposedData decomp {param.ksLevel};
    std::vector<Torus> data (10);
    initCoeffsViaUniformDistribution(data, TORUS_MIN, TORUS_MAX);
//    data = {-10,-9,9,10};
    for (auto d : data) {
        signedGadgetDecomposition(decomp, d, param);
        printf("in: %d, ", d);
        printArray(decomp.value, "decomp");
        auto recons = recomposeSelf(decomp, param);
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
    initCoeffsViaUniformDistribution(data, TORUS_MIN, TORUS_MAX);
    for (auto d : data) {
        gadgetDecompose(out, d, param);
        printf("in: %d, ", d);
        printArray(out.value, "decomp");
        auto z = recomposeSelf(out, param);
        ASSERT_EQ(z, d);
    }
    printBanner("Decompose");
}

TEST(DecompositionTest, DecomposeOverBSingleStage) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.ksLevel = 8;
    Torus mult = genIntUniformDist(TORUS_MIN, TORUS_MAX);
    std::vector<Integer> rhs(param.ksLevel);
    decomposeOverB(rhs, mult, param);
    printArray(rhs, to_string(mult) + " decomposeOverB");
    DecomposedData decomp {param.ksLevel};
    int ti = 0;
    while (ti++ < 10) {
        Torus data = genIntUniformDist(TORUS_MIN, TORUS_MAX);
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
        recompRhs[i] = recomposeSelf(decompRhs[i], param);
        ASSERT_EQ(recompRhs[i], rhs[i]);
    }
    printArray(recompRhs, "recompRhs");

    DecomposedData decompL1 {param.ksLevel};
    DecomposedData recompL1 {param.ksLevel};
//    Torus data = genIntUniformDist(TorusMin, TorusMax);
    Torus data = TORUS_MAX - 3;

//  signedGadgetDecomposition(decomp, data, param); // both correct
    gadgetDecompose(decompL1, data, param); // both correct
    printf("in: %d, ", data);
    printArray(decompL1.value, "decomp");

    // recomp first level
    recomposeFirstHalf(recompL1, decompL1, decompRhs);
    printArray(recompL1.value, "recompL1");

    // recomp second level
    auto out = recomposeSelf(recompL1, param); // equivalent to recomposeTwoParts(recompL1, decompOneOverR)
    printf("out = %d, data * mult = %d\n", out, data * mult);

    for (auto i = 0; i < param.ksLevel; i++) {
        cout << decompL1.value[i] * rhs[i] << ", "<< ((decompL1.value[i] * rhs[i]) >> (32 - 4* (i+1))) << endl;
    }

    ASSERT_EQ(out, data * mult);

    printBanner("DecomposeOverBMultiStages");
}

TEST(DecompositionTest, DecomposeTrlweTest) {
    YatfheParameters param {};
    initYatfhe(param);

    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};
    TrlweDft recompDft {param.k, param.N};

    Integer plain = 1;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    genTrlweKey(trlweKey);
//    Trlwe in2 {param.k, param.N};
//    symEncTrlweSingleSample(in2, trlweKey, mu);
    symEncTrlweSingleSampleNtt(in, inDft, trlweKey, mu);
    Trlwe intt {param.k, param.N};
    applyInttForAB(intt, inDft);
    printTrlweAB(in, "in");
//    printTrlweAB(intt, "intt");

    IntPolynomial dec2 {param.N};
    symDecTrlweToInt(dec2, in, trlweKey, param.torusBase);
    printArray(dec2.coeffs, "d2");

    // decompose
    DecomposedTrlwe out {param, 2};
//    DecomposedTrlweDft outDft {param, param.lDft}; // lDft > l, optimize?
    gadgetDecomposeTrlwe(out, in, param);
//    gadgetDecomposeTrlweNtt(outDft, inDft, param);
//    printDecomposedTrlweAB(out, "out");
//    printDecomposedTrlweNttAB(outDft, "outNtt");

    // recompose original
    recomposeTrlwe(recomp, out, param);
    printTrlweAB(in, "original");
    printTrlweAB(recomp, "recomp");
    IntPolynomial dec1 {param.N};
    symDecTrlweToInt(dec1, in, trlweKey, param.torusBase);
    for (auto j = 0; j < in.b.N; j++) {
        ASSERT_EQ(dec1.coeffs[j], plain);
    }

//    // recompose ntt
//    recomposeTrlweNtt(recompDft, outDft, param);
//    printTrlweDftAB(inDft, "inDft");
//    printTrlweDftAB(recompDft, "recompDft");
//    symDecTrlweToIntNtt(dec1, recompDft, trlweKey, param.torusBase);
//    for (auto j = 0; j < in.b.N; j++) {
//        ASSERT_EQ(dec1.coeffs[j], plain);
//    }
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
        auto z = recomposeSelf(dr, param);
        printf("a + b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a + b);
        for (auto i = 0; i < param.ksLevel; i++) {
            dr.value[i] = da.value[i] * da.sign - db.value[i] * db.sign;
        }
        z = recomposeSelf(dr, param);
        printf("a - b: a: %d, b: %d. decomp: %d, ori: %d\n", a, b, z, a - b);
    }
    printBanner("DecomposedAddSub");
}


// Actually, identical logic with DecomposeOverB test.
TEST(DecompositionTest, DecomposedMult) {
    YatfheParameters param {};
    param.q = Q_32;
    initYatfhe(param);
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
        a = genIntUniformDist(TORUS_MIN, TORUS_MAX);
        b = genIntUniformDist(TORUS_MIN, TORUS_MAX);

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
    initYatfhe(param);

    Trlwe in {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    Trlwe recomp {param.k, param.N};

    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    Torus mu = doubleToTorus32(1.0 / param.torusBase);
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


////    // decomp -> NTT -> arith -> INTT -> recomp
//    // todo
//    DecomposedTrlwe out2 {param};
//    DecomposedTrlwe out2Intt {param};
//    DecomposedTrlweDft outDft2 {param, param.l};
//    DecomposedTrlweDft multRes {param, param.l};
//    TrlweDft recompDft2 {param.k, param.N};
//    Trlwe recomp2 {param.k, param.N};
//
//    COUNT_TIME("decomp", gadgetDecomposeTrlwe(out2, in, param);) // decompose
//    COUNT_TIME("NTT",
//               for (auto i = 0; i < out2.l; i++) {
//                   applyNttForAB(outDft2.rlweDfts[i], out2.rlwes[i]); // NTT
//   })
//
//    // trgsw enc
//    Trgsw trgsw {param};
//    TrgswDft trgswDft {param};
//    Integer mu2 = 0;
//    trgswEncryptNtt(trgsw, trgswDft, param, trgswKey, mu2);
//    auto k = param.k;
//    COUNT_TIME("mult", {
//        for (size_t lvl = 0; lvl < param.l; lvl++) {
//            for (auto col = 0; col < k + 1; col++) {
//                auto &curr = (col < k) ? outDft2.rlweDfts[lvl].a[col] : outDft2.rlweDfts[lvl].b;
//                for (auto col2 = 0; col2 < k + 1; col2++) {
//                    auto &out = (col2 < k) ? multRes.rlweDfts[lvl].a[col2] : multRes.rlweDfts[lvl].b;
//                    auto &curr2 = (col2 < k) ? trgswDft.trlweDftSamples[lvl][col].a[col2]
//                                             : trgswDft.trlweDftSamples[lvl][col].b;
//                    calModularInnerProductNtt(out, curr, curr2);
//                }
//            }
//        }
//    })
//    COUNT_TIME("intt",
//               for (auto i = 0; i < out2.l; i++) {
//                   applyInttForAB(out2Intt.rlwes[i], multRes.rlweDfts[i]); // iNTT
//    })
//    COUNT_TIME("recompose", recomposeTrlwe(recomp2, out2Intt, param);)
//    TorusPolynomial dec2 {param.N};
//    symDecTrlweToInt(dec2, recomp2, trlweKey, param.torusBase);
//    printArray(dec2.coeffs, "dec2");
//    printBanner("decomp -> NTT -> recomp -> INTT");
}


// decomposed NTT does not follow the normal arithmetic rules
TEST(DecompositionTest, NttDecompArithTest) {
    YatfheParameters param {};
//    param.radixBits = 4;
//    param.lDft = 2;
//    param.k = 2;
    initYatfhe(param);

    Trlwe in {param.k, param.N};
    Trlwe in1 {param.k, param.N};
    TrlweDft inDft {param.k, param.N};
    TrlweDft inDft1 {param.k, param.N};
    Trlwe recomp {param.k, param.N};

    TrlweKey trlweKey {param.k, param.N, param.rlweStdDev};
    genTrlweKey(trlweKey);

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    Torus mu1 = doubleToTorus32(1.0 / param.torusBase);
    symEncTrlweSingleSampleNtt(in, inDft, trlweKey, mu);
    symEncTrlweSingleSampleNtt(in1, inDft1, trlweKey, mu1);


    // NTT -> decomp -> recomp -> INTT
    DecomposedTrlwe out1 {param};
    DecomposedTrlweDft outDft {param, param.lDft};
    DecomposedTrlweDft outDft1 {param, param.lDft};
    DecomposedTrlweDft outDftAdd {param, param.lDft};
    TrlweDft dft {param.k, param.N};
    TrlweDft dft1 {param.k, param.N};
    TrlweDft recompDft {param.k, param.N};
    Trlwe intt {param.k, param.N};
    COUNT_TIME("NTT", applyNttForAB(dft, in);) // NTT
    COUNT_TIME("NTT1", applyNttForAB(dft1, in1);) // NTT
    COUNT_TIME("decomp", gadgetDecomposeTrlweNtt(outDft, dft, param);) //decomp
    COUNT_TIME("decomp1", gadgetDecomposeTrlweNtt(outDft1, dft1, param);) //decomp

    TrlweDft tmp{param.k, param.N};
    DecomposedTrlweDft tmpD {param, param.lDft};
    addTrlweNtt(tmp, dft, dft1);
    gadgetDecomposeTrlweNtt(tmpD, tmp, param);

    // decomp mod add
    for (auto l = 0; l < param.lDft; l++) {
        for (auto i = 0; i < param.k; i++) {
            for (auto j = 0; j < param.N; j++) {
                outDftAdd.rlweDfts[l].a[i].coeffs[j] = (outDft.rlweDfts[l].a[i].coeffs[j] + outDft1.rlweDfts[l].a[i].coeffs[j]) % 256;
            }
        }
        for (auto j = 0; j < param.N; j++) {
            outDftAdd.rlweDfts[l].b.coeffs[j] = (outDft.rlweDfts[l].b.coeffs[j] + outDft1.rlweDfts[l].b.coeffs[j]) % 256;
        }
    }
    printArray(tmpD.rlweDfts[7].b.coeffs, "tmpD");
    printArray(outDft.rlweDfts[7].b.coeffs, "decomp");
    printArray(outDft1.rlweDfts[7].b.coeffs, "decomp1");
    printArray(outDftAdd.rlweDfts[7].b.coeffs, "outDftAdd");

    COUNT_TIME("recomp", recomposeTrlweNtt(recompDft, outDftAdd, param);) // recompose
    COUNT_TIME("INTT", applyInttForAB(intt, recompDft);) // intt

    TorusPolynomial dec1 {param.N};
    symDecTrlweToInt(dec1, intt, trlweKey, param.torusBase);
    printArray(dec1.coeffs, "dec1");
}