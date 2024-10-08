//
// Created by Xintong Song on 2024/5/28.
//
#include "yatfhe/ntt.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

void conv(const int N, const int k) {
    vector<LagrangePolynomial> a(k, LagrangePolynomial(N));
    vector<LagrangePolynomial> b(k, LagrangePolynomial(N));
    vector<Ntt14Polynomial> a14(k, Ntt14Polynomial(N));
    vector<Ntt14Polynomial> b14(k, Ntt14Polynomial(N));

    LagrangePolynomial tmpMul{N};
    Ntt14Polynomial tmpMul14{N};

    vector<TorusPolynomial> poly0(k, TorusPolynomial(N));
    vector<TorusPolynomial> poly2(k, TorusPolynomial(N));
    TorusPolynomial resMul{N};
    for (int t = 0; t < 1; ++t) {
        for (auto i = 0; i < k; i++) {
            for (auto j = 0; j < N; j++) {
                poly0[i].coeffs[j] = genIntUniformDist(CHAR_MIN, CHAR_MAX);
                poly2[i].coeffs[j] = genIntUniformDist(0, 1);
            }
        }

        COUNT_TIME("NTT_CONV64", {
            for (auto i = 0; i < k; i++) {
                applyNtt(a[i], poly0[i]);
                applyNtt(b[i], poly2[i]);
            }
            for (auto i = 0; i < k; i++) {
                calModularInnerProductNtt(tmpMul, a[i], b[i]);
            }
            applyIntt(resMul, tmpMul);
        })

        COUNT_TIME("NTT_CONV14", {
            for (auto i = 0; i < k; i++) {
                applyNtt14(a14[i], poly0[i]);
                applyNtt14(b14[i], poly2[i]);
            }
            for (auto i = 0; i < k; i++) {
                calModularInnerProductNtt14(tmpMul14, a14[i], b14[i]);
            }
            applyIntt14(resMul, tmpMul14);
        })
    }
}

int main() {
    YatfheParameters param {};
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n\n", param.n, param.k, param.N, param.radixBits, param.l);

    Ntt64 a = 1234567890123456789ULL;
    Ntt64 b = 9876543210987654321ULL;
    Ntt16 a1 = 59923;
    Ntt16 b1 = 65535;
    Ntt14 a14 = 11280;
    Ntt14 b14 = 11111;
    Ntt24 a24 = 16777216;
    Ntt24 b24 = 16777216;

    int n = param.N * param.l * param.n;
    cout << "test on " << n << " samples" << endl;
    COUNT_TIME("modMul", {
        for (auto i = 0; i < n; i++) { modMULT64(a, b); }
    })
    COUNT_TIME("modmul", {
        for (auto i = 0; i < n; i++) { modmul64(a, b); }
    })
    COUNT_TIME("fastmm", {
        for (auto i = 0; i < n; i++) { fastmm(a, b); }
    })
    COUNT_TIME("fastmm_opt", {
        for (auto i = 0; i < n; i++) { fastmm_opt(a, b); }
    })
    COUNT_TIME("16", {
        for (auto i = 0; i < n; i++) { modMULT16(a1, b1); }
    })
    COUNT_TIME("16R", {
        for (auto i = 0; i < n; i++) { modMULT16R(a1, b1); }
    })

    COUNT_TIME("14", {
        for (auto i = 0; i < n; i++) { modMULT14(a14, b14); }
    })

    COUNT_TIME("24", {
        for (auto i = 0; i < n; i++) { modMULT24(a24, b24); }
    })

    conv(param.N, param.k);

    const int N = 1024;
    IntPolynomial poly0{N};
    IntPolynomial poly2{N};;
    IntPolynomial navMul{N};
    for (auto j = 0; j < N; j++) {
        poly0.coeffs[j] = genIntUniformDist(1 << 26, 1<< 27);
        poly2.coeffs[j] = genIntUniformDist(1 << 26, 1<< 27);
    }
    COUNT_TIME("NAIVE_MULT",
               for (auto i = 0; i < 600; i++) { polynomialMulNaiveModQ(navMul, poly0, poly2, 1l<<32);})
    COUNT_TIME("NAIVE_MULT32",
               for (auto i = 0; i < 600; i++) { polynomialMulAccNaiveI32(navMul, poly0, poly2);})
}