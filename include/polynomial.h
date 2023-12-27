//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_POLYNOMIAL_H
#define HLS_YATFHE_POLYNOMIAL_H
#include "torus.h"

struct LagrangePolynomial {
    double* coeffs{ new double };
    int N{};
};

void initTorusPolynomial(TorusPolynomial& torusPolynomial, int N);

void initLagrangePolynomial(LagrangePolynomial& lagrangePolynomial, int N);

void deletePolynomial(TorusPolynomial& polynomial);

void deletePolynomial(LagrangePolynomial& polynomial);

#endif //HLS_YATFHE_POLYNOMIAL_H
