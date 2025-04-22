//
// Created by xintong on 4/21/25.
//
#include "yatfhe/blind_rotate.h"
#include "yatfhe/cmux.h"
#include <thread>

void blindRotateNormal(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = Trlwe{param.k, param.N};
        controlMux(temp, accum, input.a[i], bsk[i], param);
        accum = std::move(temp); // Update acc
    }
}

void blindRotateNormalNtt(Trlwe& accum, const vector<TrgswDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = Trlwe{param.k, param.N};
        controlMuxNtt(temp, accum, input.a[i], bskDft[i], param);
        accum = std::move(temp);
    }
}

//todo
void blindRotateGroup2(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    Trgsw tmp1{param}, tmp2{param}, tmp3{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    for (auto i = 0; i < param.n; i = i + 2) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto bsk1 = bsk[j];
        auto bsk2 = bsk[j + 1];
        auto bsk3 = bsk[j + 2];
        auto& bsk4 = bsk[j + 3];
#ifdef BLINDROT_GROUP_NAIVE
//        rotateTrgsw(bsk1, a1, param);
//        rotateTrgsw(bsk3, a1, param);
//        addTrgsw(tmp1, bsk1, bsk2);
//        addTrgsw(tmp2, bsk3, bsk4);

        std::thread t1([&](){
            rotateTrgsw(bsk1, a1, param);
            addTrgsw(tmp1, bsk1, bsk2);
        });
        std::thread t2([&](){
            rotateTrgsw(bsk3, a1, param);
            addTrgsw(tmp2, bsk3, bsk4);
        });
        t1.join();
        t2.join();

        rotateTrgsw(tmp1, a2, param);
        addTrgsw(tmp3, tmp1, tmp2);
#else
        auto a12 = a1 + a2;
        std::thread t1([&]() {
            rotateTrgsw(bsk1, a12, param);
        });
        std::thread t2([&]() {
            rotateTrgsw(bsk2, a2, param);
        });
        std::thread t3([&]() {
            rotateTrgsw(bsk3, a1, param);
        });
        t1.join();
        t2.join();
        t3.join();
//        rotateTrgsw(bsk1, a12, param);
//        rotateTrgsw(bsk2, a2, param);
//        rotateTrgsw(bsk3, a1, param);

#pragma omp parallel for collapse(3)
        for (size_t l = 0; l < param.l; l++) {
            for (size_t k = 0; k < param.k + 1; k++) {
                for (auto k2 = 0; k2 < param.k; k2++) {
#pragma omp simd
                    for (int n = 0; n < param.N; n++) {
                        tmp3.trlweSamples[l][k].a[k2].coeffs[n] = bsk1.trlweSamples[l][k].a[k2].coeffs[n] + bsk2.trlweSamples[l][k].a[k2].coeffs[n] +
                                bsk3.trlweSamples[l][k].a[k2].coeffs[n] + bsk4.trlweSamples[l][k].a[k2].coeffs[n];
                    }
                }
#pragma omp simd
                for (int n = 0; n < param.N; n++) {
                    tmp3.trlweSamples[l][k].b.coeffs[n] = bsk1.trlweSamples[l][k].b.coeffs[n] + bsk2.trlweSamples[l][k].b.coeffs[n] +
                                                              bsk3.trlweSamples[l][k].b.coeffs[n] + bsk4.trlweSamples[l][k].b.coeffs[n];
                }
            }
        }
#endif
        temp = Trlwe{param.k, param.N};
        externalProductTrgsw(temp, tmp3, accum, param);
        accum = std::move(temp);
        j += batchSize;
    }
}

void blindRotateGroup2Ntt(Trlwe& accum, const vector<TrgswDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    TrgswDft tmp1{param}, tmp2{param}, tmp3{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    for (auto i = 0; i < param.n; i = i + 2) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto bsk1 = bskDft[j];
        auto& bsk2 = bskDft[j + 1];
        auto bsk3 = bskDft[j + 2];
        auto& bsk4 = bskDft[j + 3];

        rotateTrgswNtt(bsk1, a1, param);
        rotateTrgswNtt(bsk3, a1, param);
        addTrgswNtt(tmp1, bsk1, bsk2);
        addTrgswNtt(tmp2, bsk3, bsk4);

        rotateTrgswNtt(tmp1, a2, param);
        addTrgswNtt(tmp3, tmp1, tmp2);

        temp = Trlwe{param.k, param.N};
        externalProductTrgswNtt(temp, tmp3, accum, param);
        accum = std::move(temp);
        j += batchSize;
    }
}

/**
 * Multiply the accumulator by X^sum(bara_i * s_i)
 * */
void blindRotate(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    switch(param.group) {
        case 2:
            blindRotateGroup2(accum, bsk, input, param);
            return;
        default:
            blindRotateNormal(accum, bsk, input, param);
    }
}

void blindRotateNtt(Trlwe& accum, const vector<TrgswDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    switch(param.group) {
        case 2:
            blindRotateGroup2Ntt(accum, bskDft, input, param);
            return;
        default:
            blindRotateNormalNtt(accum, bskDft, input, param);
    }
}

void blindRotateApproxCRT(std::vector<Trlwe8>& accum, const vector<vector<Trgsw8>>& bskCRT, const ScaledTlwe& input, const YatfheParameters& param) {
    std::vector<Trlwe8> temp(param.d, Trlwe8{param.k, param.N});
    for (size_t i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        for (size_t d = 0; d < param.d; d++) {
            temp[d] = Trlwe8{param.k, param.N};
        }
        controlMuxApproxCRT(temp, accum, input.a[i], bskCRT[i], param);
        swap(accum, temp);
    }
}

void blindRotateApproxCRTNtt(std::vector<Trlwe8>& accum, const vector<vector<TrgswDft24>>& bskCRT, const ScaledTlwe& input, const YatfheParameters& param) {
    std::vector<Trlwe8> temp(param.d, Trlwe8{param.k, param.N});
    for (size_t i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        for (size_t d = 0; d < param.d; d++) {
            temp[d] = Trlwe8{param.k, param.N};
        }
        controlMuxApproxCRTNtt(temp, accum, input.a[i], bskCRT[i], param);
        swap(accum, temp);
    }
}
