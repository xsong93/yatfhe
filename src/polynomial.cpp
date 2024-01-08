//
// Created by Xintong Song on 2023/12/25.
//
#include "polynomial.h"

void initTorusPolynomial(TorusPolynomial& torusPolynomial, const int N) {
    torusPolynomial.N = N;
}

void initLagrangePolynomial(LagrangePolynomial& lagrangePolynomial, const int N) {
    lagrangePolynomial.N = N;
}

void deletePolynomial(TorusPolynomial& polynomial) {
//    delete polynomial.coeffs;
//    polynomial.coeffs = nullptr;
}

void deletePolynomial(LagrangePolynomial& polynomial) {
    delete polynomial.coeffs;
    polynomial.coeffs = nullptr;
}