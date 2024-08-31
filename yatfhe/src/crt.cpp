//
// Created by Xintong Song on 2024/8/27.
//
#include "yatfhe/crt.h"
#include "yatfhe/numeric_functions.h"

int calApproxCRTError(const std::vector<int>& f_tilde, const std::vector<int>& coeffs) {
    int infNorm = 0;
    for (auto i = 0; i < f_tilde.size(); i++) {
        infNorm  = std::max(infNorm, std::abs(f_tilde[i] - coeffs[i]));
    }
    return infNorm;
}

void calGadgetVector(std::vector<long>& w, const int Qlow, const std::vector<int>& highModuli) {
    int Q = highModuli[0];
    for (auto i = 1; i < highModuli.size(); i++) {
        Q *= highModuli[i];
    }
    auto Q1 = Q / highModuli[0];
    auto Q2 = Q / highModuli[1];
    w[0] = Qlow * Q1 * (modInverse(Qlow * Q1, highModuli[0]));
    w[1] = Qlow * Q2 * (modInverse(Qlow * Q2, highModuli[1]));
//    printf("w1:%ld, w2:%ld\n", w[0], w[1]);
}

void approximateCRTDecomposition(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const int Qlow, const std::vector<int>& lowModuli, const std::vector<int>& highModuli) {
    for (int i = 0; i < coeffs.size(); ++i) {
        int fi = coeffs[i];

        // Calculate the low part contribution
        int lowSum = 0;
        for (int u = 0; u < lowModuli.size(); ++u) {
            int q_u = lowModuli[u];
            int inv = (int)modInverse(Qlow / q_u, q_u);
            lowSum += (int)intModP(Qlow / q_u * intModP(inv * intModP(fi, q_u), q_u), Qlow);
        }

        // Calculate f for high moduli
        for (auto j = 0; j < f.size(); j++) {
            f[j][i] = (int8_t) intModP(fi - lowSum, highModuli[j]);
        }
    }
}

void approximatePolynomialReconstruction(std::vector<int>& f_tilde, const std::vector<std::vector<int8_t>>& f, const std::vector<long>& w, long q) {
    for (int i = 0; i < f_tilde.size(); i++) {
        for (int j = 0; j < w.size(); j++) {
            f_tilde[i] = (int) intModP(intModP((long)f[j][i] * w[j], q) + f_tilde[i], q);
        }
    }
}