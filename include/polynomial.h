//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H
#include "torus.h"

struct LagrangePolynomial {
    double* coeffs;
    int N;
};

//TorusPolynomial* initTorusPolynomial(int N);void initTorusPolynomial(TorusPolynomial* torusPolynomial, int N)

//LagrangePolynomial* initLagrangePolynomial(int N);

void initTorusPolynomial(TorusPolynomial* torusPolynomial, int N);

void initLagrangePolynomial(LagrangePolynomial* lagrangePolynomial, int N);


#endif //HLS_YATFHE_POLYNOMIAL_H
