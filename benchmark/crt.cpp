//
// Created by Xintong Song on 2024/8/27.
//
#include "yatfhe/crt.h"
#include "yautil/tool.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/initializer.h"

int main() {
    YatfheParameters p {};
    yatfheInit(p);
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

    printArray(p.taoU, "taoU");
    printArray(p.taoUInv, "taoUInv");
    cout << longModP(p.taoU[2] * p.taoUInv[2], p.qd[2]) << endl;
    cout << longModP(p.taoU[3] * p.taoUInv[3], p.qd[3]) << endl;

/*    std::vector<std::vector<std::vector<long>>> iii(1000, std::vector<std::vector<long>>(1000, std::vector<long>(1000, 0)));

    COUNT_TIME("n3",
               for (size_t i = 0; i < 1000; i++) {
                   auto& t11 = iii[i];
                   for (size_t j = 0; j < 1000; j++) {
                       auto& t12 = t11[j];
                       for (size_t k = 0; k < 1000; k++) {
                           t12[k] = i * j * k;
                       }
                   }
               })
    COUNT_TIME("n21",
               for (size_t i = 0; i < 1000; i++) {
                   auto& t11 = iii[i];
                   for (size_t j = 0; j < 1000; j++) {
                       for (size_t k = 0; k < 1000; k++) {
                           t11[j][k] = i * j * k;
                       }
                   }
               })
    COUNT_TIME("n22",
               for (size_t i = 0; i < 1000; i++) {
                   for (size_t j = 0; j < 1000; j++) {
                       auto& t12 = iii[i][j];
                       for (size_t k = 0; k < 1000; k++) {
                           t12[k] = i * j * k;
                       }
                   }
               })
    COUNT_TIME("n1",
               for (size_t i = 0; i < 1000; i++) {
                   for (size_t j = 0; j < 1000; j++) {
                       for (size_t k = 0; k < 1000; k++) {
                           iii[i][j][k] = i * j * k;
                       }
                   }
               })*/


}