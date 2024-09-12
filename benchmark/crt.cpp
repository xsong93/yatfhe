//
// Created by Xintong Song on 2024/8/27.
//
#include "yatfhe/crt.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"

int main() {
    std::vector<int> coeffs = {656381177, -1322693974, 749894848, 1618033988};
    int Qlow = 55687;
    int l = 2;
    std::vector<int> lowModuli = {233, 239};
    std::vector<int> highModuli = {241, 251};
    long q = 3368562317;
    std::vector<std::vector<int8_t>> f(l, std::vector<int8_t>(coeffs.size(), 0));
    std::vector<int> f_tilde(coeffs.size(), 0);
    std::vector<long> w(l, 0l);

    approxCRTDecomp(f, coeffs, Qlow, lowModuli, highModuli);
    calGadgetVector(w, Qlow, highModuli);
    approxPolyReconstruct(f_tilde, f, w, q);
    for (auto i = 0; i < f.size(); i++) {
        printArray(f[i], "f" + to_string(i) + "(mod " + to_string(highModuli[i]) + ")");
    }
    printArray(coeffs, "f_origi");
    printArray(f_tilde, "f_tilde");
    cout << "||f - f_t|| = " << calApproxCRTError(f_tilde, coeffs) << endl;

    YatfheParameters p {};
    vector<Trlwe8> t1 (p.d, Trlwe8(p.k, p.N));
    vector<Trlwe8> t2 (p.dh, Trlwe8(p.k, p.N));

    for (size_t i = 0; i < p.d; i++) {
        for (size_t j = 0; j < p.N; j++) {
            for (size_t k = 0; k < p.k; k++) {
                t1[i].a[k].coeffs[j] = (int8_t)intModP(j, 256);
            }
            t1[i].b.coeffs[j] = (int8_t)intModP(j, 256);
        }
    }
    syncGadgetDecomp(t2, t1, p);
//    for (size_t i = 0; i < p.dh; i++) {
//        for (size_t j = 0; j < p.N; j++) {
//            for (size_t k = 0; k < p.k; k++) {
//                printf("t2a: %d; ", t2[i].a[k].coeffs[j]);
//            }
//            cout<< endl << endl;
//            printf("t2b: %d; ", t2[i].b.coeffs[j]);
//        }
//    }
    cout << longModP(18014398509481984, 1l << 32);
}