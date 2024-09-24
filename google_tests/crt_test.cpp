//
// Created by Xintong Song on 2024/9/24.
//
#include "gtest/gtest.h"
#include "yatfhe/crt.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/initializer.h"

TEST(CRT, APPROX_CRT) {
    YatfheParameters p {};
    yatfheInit(p);
    std::vector<int> coeffs = {656381177, -1322693974, 749894848, 1618033988};
    int Qlow = 55687;
    int l = 2;
    std::vector<int> lowModuli = {239, 233};
    std::vector<int> highModuli = {251, 241};
    long q = 3368562317;
    std::vector<std::vector<int8_t>> f(l, std::vector<int8_t>(coeffs.size(), 0));
    std::vector<int> f_tilde(coeffs.size(), 0);
    std::vector<long> w(l, 0l);

    approxCRTDecomp(f, coeffs, Qlow, lowModuli, highModuli);
    calGadgetVector(w, Qlow, highModuli);
    approxCRTReconstructPoly(f_tilde, f, w, q);
    for (auto i = 0; i < f.size(); i++) {
        printArray(f[i], "f" + to_string(i) + "(mod " + to_string(highModuli[i]) + ")");
    }
    printArray(coeffs, "f_origi");
    printArray(f_tilde, "f_tilde");
    int err = calApproxCRTError(f_tilde, coeffs);
    cout << "||f - f_t|| = " << err << endl;
    int maxErr = p.dl * p.qLow / 2;
    ASSERT_LE(err, maxErr);
    printBanner("APPROX_CRT");
}

TEST(CRT, EXACT_CRT) {
    YatfheParameters p {};
    yatfheInit(p);
    std::vector<int32_t> coeffs = {656381177, -1322693974, 749894848, 1618033988};
    std::vector<std::vector<int8_t>> f(p.d, std::vector<int8_t>(coeffs.size(), 0));
    std::vector<int32_t> f_tilde(coeffs.size(), 0);
    exactCRTDecomp(f, coeffs, p);
    for (auto i = 0; i < f.size(); i++) {
        printArray(f[i], "f" + to_string(i) + "(mod " + to_string(p.qd[i]) + ")");
    }
    exactCRTReconstruct(f_tilde, f, p);
    printArray(coeffs, "f_origi");
    printArray(f_tilde, "f_tilde");
    for (size_t i = 0; i < coeffs.size(); i++) {
        ASSERT_EQ(coeffs[i], f_tilde[i]);
    }
    printBanner("EXACT_CRT");
}