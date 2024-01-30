//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "ntt.h"
#include "polynomial.h"
#include "numeric_functions.h"

TEST(NttTest, X) {
    LagrangePolynomial a(1024);
    LagrangePolynomial b(1024);
    LagrangePolynomial res(1024);
    LagrangePolynomial res1(1024);
    TorusPolynomial a1(1024);
    TorusPolynomial b1(1024);
    TorusPolynomial res2(1024);
    for (int i = 0; i < a.N; i++) {
        a1.coeffs[i] = i;
        b1.coeffs[i] = i;
    }
    applyNtt(a, a1);
    applyNtt(b, b1);
    for (int i = 0; i < a.N; i++) {
        res.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(res1, res);
    uint64_t med = MODULUS / 2;
    std::cout << "res1:[";
    for (int i = 0; i < res2.N; i++) {
        res2.coeffs[i] = (Torus)((res1.coeffs[i] & 0xffffffff) - (res1.coeffs[i] > med));
        std::cout << i << ":" <<res2.coeffs[i]<<" ";
    }
    std::cout<<std::endl;
}