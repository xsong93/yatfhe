//
// Created by Xintong Song on 2023/12/25.
//
#include "polynomial.h"

//TorusPolynomial* initTorusPolynomial(int N) {
//    TorusPolynomial* torusPolynomial{new TorusPolynomial};
//    torusPolynomial->coeffs = new Torus;
//    torusPolynomial->N = N;
//    return torusPolynomial;
//}

void initTorusPolynomial(TorusPolynomial* torusPolynomial, const int N) {
    torusPolynomial->coeffs = new Torus;
    torusPolynomial->N = N;
}

//LagrangePolynomial* initLagrangePolynomial(int N) {
//    LagrangePolynomial* lagrangePolynomial {new LagrangePolynomial};
//    lagrangePolynomial->coeffs = new double;
//    lagrangePolynomial->N = N;
//    return lagrangePolynomial;
//}

void initLagrangePolynomial(LagrangePolynomial* lagrangePolynomial, const int N) {
    lagrangePolynomial->coeffs = new double;
    lagrangePolynomial->N = N;
}