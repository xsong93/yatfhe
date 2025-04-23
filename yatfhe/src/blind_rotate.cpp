//
// Created by xintong on 4/21/25.
//
#include "yatfhe/blind_rotate.h"
#include "yatfhe/cmux.h"
#include "yatfhe/ntt_hexl.h"
#include "yautil/time_counter.h"
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

        for (size_t l = 0; l < param.l; l++) {
            for (size_t k = 0; k < param.k + 1; k++) {
                for (auto k2 = 0; k2 < param.k; k2++) {
                    auto& coeffA1 = bsk1.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA2 = bsk2.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA3 = bsk3.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA4 = bsk4.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffAT = tmp3.trlweSamples[l][k].a[k2].coeffs;
                    for (int n = 0; n < param.N; n++) {
//                        tmp3.trlweSamples[l][k].a[k2].coeffs[n] = bsk1.trlweSamples[l][k].a[k2].coeffs[n]
//                                + bsk2.trlweSamples[l][k].a[k2].coeffs[n]
//                                + bsk3.trlweSamples[l][k].a[k2].coeffs[n]
//                                + bsk4.trlweSamples[l][k].a[k2].coeffs[n];
                        coeffAT[n] = coeffA1[n] + coeffA2[n] + coeffA3[n] + coeffA4[n];
                    }
                }
                auto& coeffB1 = bsk1.trlweSamples[l][k].b.coeffs;
                auto& coeffB2 = bsk2.trlweSamples[l][k].b.coeffs;
                auto& coeffB3 = bsk3.trlweSamples[l][k].b.coeffs;
                auto& coeffB4 = bsk4.trlweSamples[l][k].b.coeffs;
                auto& coeffBT = tmp3.trlweSamples[l][k].b.coeffs;
                for (int n = 0; n < param.N; n++) {
//                    tmp3.trlweSamples[l][k].b.coeffs[n] = bsk1.trlweSamples[l][k].b.coeffs[n]
//                            + bsk2.trlweSamples[l][k].b.coeffs[n]
//                            + bsk3.trlweSamples[l][k].b.coeffs[n]
//                            + bsk4.trlweSamples[l][k].b.coeffs[n];
                    coeffBT[n] = coeffB1[n] + coeffB2[n] + coeffB3[n] + coeffB4[n];
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
        auto bsk2 = bskDft[j + 1];
        auto bsk3 = bskDft[j + 2];
        auto& bsk4 = bskDft[j + 3];
#ifdef BLINDROT_GROUP_NAIVE
        std::thread t1([&](){
            rotateTrgswNtt(bsk1, a1, param);
            addTrgswNtt(tmp1, bsk1, bsk2);
        });
        std::thread t2([&](){
            rotateTrgswNtt(bsk3, a1, param);
            addTrgswNtt(tmp2, bsk3, bsk4);
        });
        t1.join();
        t2.join();

        rotateTrgswNtt(tmp1, a2, param);
        addTrgswNtt(tmp3, tmp1, tmp2);
#else
        auto a12 = a1 + a2;
        std::thread t1([&]{
            rotateTrgswNtt(bsk1, a12, param);
        });
        std::thread t2([&]{
            rotateTrgswNtt(bsk2, a2, param);
        });
        std::thread t3([&]{
            rotateTrgswNtt(bsk3, a1, param);
        });
        t1.join();
        t2.join();
        t3.join();
        auto q = NttHexl::getNttHexl().GetModulus();
        auto twiceMod = 2 * q;
        auto L = param.l;
        auto K = param.k;
        auto N = param.N;
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K + 1; k++) {
                for (int k2 = 0; k2 < K; k2++) {
                    auto& coeffA1 = bsk1.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA2 = bsk2.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA3 = bsk3.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA4 = bsk4.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffAT = tmp3.trlweDftSamples[l][k].a[k2].coeffs;
                    for (int n = 0; n < N; n++) {
//                        auto tmp = coeffA1[n] + coeffA2[n] + coeffA3[n] + coeffA4[n];
//                        coeffAT[n] = ReduceMod<4>(tmp, q, &twiceMod);
                        auto val1 = AddUIntMod(coeffA1[n], coeffA2[n], q);
                        auto val2 = AddUIntMod(coeffA3[n], coeffA4[n], q);
                        coeffAT[n] = AddUIntMod(val1, val2, q);
                    }
                }
                auto& coeffB1 = bsk1.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB2 = bsk2.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB3 = bsk3.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB4 = bsk4.trlweDftSamples[l][k].b.coeffs;
                auto& coeffBT = tmp3.trlweDftSamples[l][k].b.coeffs;
                for (int n = 0; n < param.N; n++) {
//                    auto tmp = coeffB1[n] + coeffB2[n] + coeffB3[n] + coeffB4[n];
//                    coeffBT[n] = ReduceMod<4>(tmp, q, &twiceMod);
                    uint64_t val1 = AddUIntMod(coeffB1[n], coeffB2[n], q);
                    uint64_t val2 = AddUIntMod(coeffB3[n], coeffB4[n], q);
                    coeffBT[n] = AddUIntMod(val1, val2, q);
                }
            }
        }
#endif
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
