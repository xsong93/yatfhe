//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

using namespace NttNative;
using namespace NttNative64;

std::vector<int> extractValues(const std::string& input) {
    std::vector<int> result;
    std::istringstream iss(input);
    std::string token;

    while (std::getline(iss, token, ' ')) {
        size_t pos = token.find(':');
        if (pos != std::string::npos) {
            std::string valueStr = token.substr(pos + 1);
            int value = std::stoi(valueStr);
            result.push_back(value);
        }
    }

    return result;
}

TEST(NTT64, NTT_INTT) {
    YatfheParameters p {};
    p.N = 512;
    initYatfhe(p);

    auto N = p.N;
    NttPolynomial resNtt{N};
    TorusPolynomial a1{N};
//    TorusPolynomial resIntt{N};
    int t = 10;
    while (t-- > 0) {
        for (auto i = 0; i < N; i++) {
            a1.coeffs[i] = genIntUniformDist(INT_MIN_VALUE, INT_MAX_VALUE);
        }
        TorusPolynomial resIntt{a1};
        int t2 = 100;
        while (t2 -- > 0) {
            applyNtt(resNtt, resIntt);
            applyIntt(resIntt, resNtt);
        }
        printArray(resNtt.coeffs, "resNtt");
        printArray(a1.coeffs, "orig");
        printArray(resIntt.coeffs, "intt");
        for (auto i = 0; i < a1.N; i++) {
            if (a1.coeffs[i] != resIntt.coeffs[i]) {
                printf("NE at index %d\n", i);
            }
            ASSERT_EQ(a1.coeffs[i], resIntt.coeffs[i]);
        }
    }
    printBanner("NTT64.NTT_INTT");
}

TEST(NTT64, ADD_CONST) {
    YatfheParameters p {};
    p.N = 1024;
    initYatfhe(p);

    auto N = p.N;
    NttPolynomial a{N};
    NttPolynomial b{N};
    NttPolynomial resNtt{N};
    IntPolynomial poly{N};
    IntPolynomial c{N};
    IntPolynomial res{N};
    for (int i = 0; i < a.N; i++) {
        poly.coeffs[i] = i;
    }
    c.coeffs[0] = 77;
    applyNtt(a, poly);
    applyNtt(b, c);
    printArray(c.coeffs, "cOri");
    printArray(b.coeffs, "bNtt");
    for (int i = 0; i < a.N; i++) {
        resNtt.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(res, resNtt);
    printArray(res.coeffs, "res");
}

TEST(NTT64, BASIC_ARITH) {
    YatfheParameters p {};
    p.N = 1024;
    // Pinned to a 32-bit torus
    p.q = Q_32;
    p.torusBits = 32;
    initYatfhe(p);

    auto N = p.N;
    NttPolynomial a{N};
    NttPolynomial b{N};
    NttPolynomial c{N};
    NttPolynomial tmpMul{N};
    NttPolynomial tmpMul1{N};
    NttPolynomial tmpAdd{N};
    NttPolynomial tmpSub{N};

    IntPolynomial poly0{N};
    IntPolynomial poly1{N};
    IntPolynomial poly2{N};
    IntPolynomial resMul{N};
    IntPolynomial resMul1{N};
    IntPolynomial resMul2{N};
    IntPolynomial resAdd{N};
    IntPolynomial resSub{N};
    IntPolynomial navMul{N};
    IntPolynomial navAdd{N};
    IntPolynomial navSub{N};
    int t = 1;
    while (t-- > 0) {
        for (auto j = 0; j < N; j++) {
            poly0.coeffs[j] = genIntUniformDist(1 << 12, 1 << 15);
            poly1.coeffs[j] = genIntUniformDist(1 << 12, 1 << 15);
            poly2.coeffs[j] = genIntUniformDist(1 << 12, 1 << 24);
        }
        printArray(poly0.coeffs, "poly0");
        printArray(poly1.coeffs, "poly1");
        printArray(poly2.coeffs, "poly2");

        COUNT_TIME("NTT_MULT_CHAINING", {
            applyNtt(a, poly0);
            applyNtt(b, poly1);
            applyNtt(c, poly2);
            for (int i = 0; i < a.N; i++) {
                auto tmp = fastmm_opt(a.coeffs[i], b.coeffs[i]);
                tmpMul.coeffs[i] = fastmm_opt(tmp, c.coeffs[i]);
            }
            applyIntt(resMul, tmpMul);
        })
        COUNT_TIME("NTT_MULT_NORMAL", {
            applyNtt(a, poly0);
            applyNtt(b, poly1);
            applyNtt(c, poly2);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = fastmm_opt(a.coeffs[i], b.coeffs[i]);
            }
            applyIntt(resMul1, tmpMul);
            applyNtt(tmpMul1 ,resMul1);
            for (int i = 0; i < a.N; i++) {
                tmpMul.coeffs[i] = fastmm_opt(tmpMul1.coeffs[i], c.coeffs[i]);
            }
            applyIntt(resMul2, tmpMul);
        })
        COUNT_TIME("NAIVE_MULT",
            IntPolynomial tmp{N};
            multIntPolynomialModQ(tmp, poly0, poly1, TORUS_Q);
            multIntPolynomialModQ(navMul, tmp, poly2, TORUS_Q);
        )

        for (int i = 0; i < a.N; i++) {
            tmpAdd.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
            tmpSub.coeffs[i] = modSub(a.coeffs[i], b.coeffs[i]);
        }
        applyIntt(resAdd, tmpAdd);
        applyIntt(resSub, tmpSub);
        addIntPolynomial(navAdd, poly0, poly1);
        subIntPolynomial(navSub, poly0, poly1);

        printArray(resMul.coeffs, "resMul");
        printArray(resMul2.coeffs, "resMul2");
        printArray(navMul.coeffs, "TRUE");

        EXPECT_EQ(resMul2.coeffs, navMul.coeffs);
        EXPECT_NE(resMul.coeffs, navMul.coeffs); // overflow
//            EXPECT_EQ(resAdd.coeffs[i], navAdd.coeffs[i]);
//            EXPECT_EQ(resSub.coeffs[i], navSub.coeffs[i]);

    }
    printBanner("NTT64.BASIC_ARITH");
}

TEST(NTT64, CONVOLUTION) {
    YatfheParameters p {};
    p.N = 1024;
    // Pinned to a 32-bit torus
    p.q = Q_32;
    p.torusBits = 32;
    initYatfhe(p);

    auto N = p.N;
    auto k = p.k;

    vector<NttPolynomial> a(k, NttPolynomial(N));
    vector<NttPolynomial> b(k, NttPolynomial(N));
    NttPolynomial tmpMul{N};

    vector<IntPolynomial> poly0(k, TorusPolynomial(N));
    vector<IntPolynomial> poly2(k, TorusPolynomial(N));
    IntPolynomial resMul{N};
    IntPolynomial navMul{N};
    for (int t = 0; t < 10; ++t) {
        for (auto i = 0 ; i < k; i++) {
            for (auto j = 0; j < N; j++) {
                poly0[i].coeffs[j] = genIntUniformDist(TORUS_MIN, TORUS_MAX);
                poly2[i].coeffs[j] = genIntUniformDist(0, 1);
            }
            printArray(poly0[i].coeffs, "poly0" + to_string(i));
            printArray(poly2[i].coeffs, "poly2" + to_string(i));
        }

        COUNT_TIME("NTT_MULT", {
            for (auto i = 0 ; i < k; i++) {
                applyNtt(a[i], poly0[i]);
                applyNtt(b[i], poly2[i]);
            }

            calModularInnerProductNtt(tmpMul, a, b);
//            printArray(tmpMul.coeffs,"tmpMUL");
            applyIntt(resMul, tmpMul);})
        COUNT_TIME("NAIVE_MULT",
                   for (auto i = 0 ; i < k; i++) {
                       multIntPolynomialAcc(navMul, poly0[i], poly2[i]);
                   })
        for (int idx = 0; idx < navMul.N; idx++) {
            navMul.coeffs[idx] = static_cast<Torus>(longModP(navMul.coeffs[idx], TORUS_Q));
        }
        printArray(resMul.coeffs, "resMul");
        printArray(navMul.coeffs, "navMul");

        for (int i = 0; i < navMul.N; i++) {
            if (resMul.coeffs[i] != navMul.coeffs[i]){
                cout<<"error at:"<<i<<endl;
            }
            EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        }
    }

    printBanner("NTT64.CONVOLUTION");
}