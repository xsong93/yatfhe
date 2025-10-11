//
// Created by xintong on 4/21/25.
//
#include <omp.h>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/cmux.h"
#include "yatfhe/ntt_hexl.h"
#include "yautil/time_counter.h"
#include "yautil/multi_threading.h"

void preRotateBinary(Trlwe& trlweOut, vector<TrgswMPDft>& trgswDftsOut, const vector<Trlwe>& key1,
                     const vector<vector<TrgswMPDft>>& trgswDftsIn, const ScaledTlwe& in, const int batchSize,
                     const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();

    vector<future<void>> futures;
    futures.reserve(batchSize);

    {
        Trlwe tmp{param.k, param.N};
        rotateTrlwe(tmp, key1[0], in.a[0]);
        addTrlwe(tmp, tmp, key1[1]);
        rotateTrlwe(trlweOut, tmp, -in.b);
    }

    for (int start = 1; start < n; start += batchSize) {
        futures.clear();
        const int end = min(start + batchSize, n);

        for (int i = start; i < end; ++i) {
            const auto ai = in.a[i];
            const int j = i - 1;
            const auto& keyIn = trgswDftsIn[j];
            auto& keyOut = trgswDftsOut[j];
            futures.emplace_back(pool.enqueue([ai, &param, &keyOut, &keyIn] {
                keyOut = keyIn[0];
                if (ai != 0) {
                    rotateTrgswMPNtt(keyOut, ai, param);
                }
                addTrgswMPNtt(keyOut, keyOut, keyIn[1]);
            }));
        }

        for (auto& f : futures) {
            f.get();
        }
    }
}

void preRotateTernary(Trlwe& trlweOut, vector<TrgswMPDft>& trgswDftsOut, const vector<Trlwe>& key1,
                      const vector<vector<TrgswMPDft>>& trgswDftsIn, const ScaledTlwe& in, const int batchSize,
                      const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();

    vector<future<void>> futures;
    futures.reserve(batchSize);

    {
        Trlwe tmp{param}, tmp2{param};
        rotateTrlwe(tmp, key1[0], in.a[0]);
        rotateTrlwe(tmp2, key1[1], -in.a[0]);
        addTrlwe(tmp, tmp, tmp2, key1[2]);
        rotateTrlwe(trlweOut, tmp, -in.b);
    }

    for (int start = 1; start < n; start += batchSize) {
        futures.clear();
        const int end = min(start + batchSize, n);

        for (int i = start; i < end; ++i) {
            const auto ai = in.a[i];
            const int j = i - 1;
            const auto& keyIn = trgswDftsIn[j];
            auto& keyOut = trgswDftsOut[j];
            futures.emplace_back(pool.enqueue([ai, &param, &keyOut, &keyIn] {
                keyOut = keyIn[0];
                auto keyOut2 = keyIn[1];
                if (ai != 0) {
                    rotateTrgswMPNtt(keyOut, ai, param);
                    rotateTrgswMPNtt(keyOut2, -ai, param);
                }
                addTrgswMPNtt(keyOut, keyOut, keyOut2, keyIn[2]);
            }));
        }

        for (auto& f : futures) {
            f.get();
        }
    }
}

void preRotateTernary2(Trlwe& trlweOut, vector<TrgswMPDft>& trgswDftsOut, const vector<Trlwe>& key1,
                      const vector<vector<TrgswMPDft>>& trgswDftsIn, const ScaledTlwe& in, const int batchSize,
                       const int tasksPerThread, const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();

    vector<future<void>> futures;
    futures.reserve(batchSize);

    {
        Trlwe tmp{param}, tmp2{param};
        rotateTrlwe(tmp, key1[0], in.a[0]);
        rotateTrlwe(tmp2, key1[1], -in.a[0]);
        addTrlwe(tmp, tmp, tmp2, key1[2]);
        rotateTrlwe(trlweOut, tmp, -in.b);
    }

    for (int start = 1; start < n; start += batchSize * tasksPerThread) {
        futures.clear();
        const int end = min(start + batchSize * tasksPerThread, n);

        // Process tasks in chunks of 'tasksPerThread'
        for (int chunkStart = start; chunkStart < end; chunkStart += tasksPerThread) {
            const int chunkEnd = min(chunkStart + tasksPerThread, end);

            futures.emplace_back(pool.enqueue([chunkStart, chunkEnd, &in, &param, &trgswDftsOut, &trgswDftsIn] {
                for (int i = chunkStart; i < chunkEnd; ++i) {
                    const auto ai = in.a[i];
                    const int j = i - 1;
                    auto& keyOut = trgswDftsOut[j];
                    const auto& keyIn = trgswDftsIn[j];

                    keyOut = keyIn[0];
                    auto keyOut2 = keyIn[1];
                    if (ai != 0) {
                        rotateTrgswMPNtt(keyOut, ai, param);
                        rotateTrgswMPNtt(keyOut2, -ai, param);
                    }
                    addTrgswMPNtt(keyOut, keyOut, keyOut2, keyIn[2]);
                }
            }));
        }

        for (auto& f : futures) {
            f.get();
        }
    }
}

void preRotateInternal(vector<TrgswMP>& trgswsOut, vector<TrgswMPDft>& trgswDftsOut, vector<vector<TrgswMP>>& trgswsIn,
    vector<vector<TrgswMPDft>>& trgswDftsIn, const std::vector<int>& aV, const YatfheParameters& param) {
    auto const n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n/2; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n/2); i++) {
            futures.emplace_back(pool.enqueue([i, &trgswsIn, &trgswDftsIn, &aV, &param, &trgswsOut, &trgswDftsOut] {
                const int j = i * 2;
                auto& rotated1 = trgswsIn[i][0];
                auto& rotated2 = trgswDftsIn[i][0];

                rotateTrgswMP(rotated1, aV[j], param);
                rotateTrgswMPNtt(rotated2, aV[j+1], param);

                addTrgswMP(trgswsOut[i], rotated1, trgswsIn[i][1]);
                addTrgswMPNtt(trgswDftsOut[i], rotated2, trgswDftsIn[i][1]);
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }

    // for (int i = 0; i < n / 2; i++) {
    //     const int j = i * 2;
    //     auto& rotated1 = trgswsIn[i][0];
    //     auto& rotated2 = trgswDftsIn[i][0];
    //
    //     rotateTrgswMP(rotated1, aV[j], param);
    //     rotateTrgswMPNtt(rotated2, aV[j+1], param);
    //
    //     addTrgswMP(trgswsOut[i], rotated1, trgswsIn[i][1]);
    //     addTrgswMPNtt(trgswDftsOut[i], rotated2, trgswDftsIn[i][1]);
    // }
}

void preRotateInternalAsym(vector<Trlev>& trglevsOut, vector<TrgswMPDft>& trgswDftsOut, vector<vector<Trlev>>& trlevsIn,
    vector<vector<TrgswMPDft>>& trgswDftsIn, const std::vector<int>& aV, const YatfheParameters& param) {
    auto const n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n/2; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n/2); i++) {
            futures.emplace_back(pool.enqueue([i, &trlevsIn, &trgswDftsIn, &aV, &param, &trglevsOut, &trgswDftsOut] {
                const int j = i * 2;
                rotateTrgswMPNtt(trgswDftsIn[i][0], aV[j], param);
                rotateTrlev(trlevsIn[i][0], aV[j + 1], param);
                addTrgswMPNtt(trgswDftsOut[i], trgswDftsIn[i][0], trgswDftsIn[i][1]);
                addTrlev(trglevsOut[i], trlevsIn[i][0], trlevsIn[i][1]);
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

void preRotateInternalAsymOpt(vector<Trlev>& trglevsOut, Trlwe& trlweOut, vector<TrgswMPDft>& trgswDftsOut,
                              vector<vector<Trlev>>& trlevsIn, vector<Trlwe>& trlwesIn, vector<vector<TrgswMPDft>>& trgswDftsIn,
                              const std::vector<int>& aV, const YatfheParameters& param) {
    auto const n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n/2; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n/2); i++) {
            if (i == n/2 - 1) {
                const int j = i * 2;
                Trlwe tmp {param.k, param.N};
                rotateTrgswMPNtt(trgswDftsIn[i][0], aV[j], param);
                rotateTrlwe(tmp, trlwesIn[0], aV[j + 1]);
                addTrgswMPNtt(trgswDftsOut[i], trgswDftsIn[i][0], trgswDftsIn[i][1]);
                addTrlwe(trlweOut, tmp, trlwesIn[1]);
                continue;
            }
            futures.emplace_back(pool.enqueue([i, &trlevsIn, &trgswDftsIn, &aV, &param, &trglevsOut, &trgswDftsOut] {
                const int j = i * 2;
                rotateTrgswMPNtt(trgswDftsIn[i][0], aV[j], param);
                rotateTrlev(trlevsIn[i][0], aV[j + 1], param);
                addTrgswMPNtt(trgswDftsOut[i], trgswDftsIn[i][0], trgswDftsIn[i][1]);
                addTrlev(trglevsOut[i], trlevsIn[i][0], trlevsIn[i][1]);
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

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

void blindRotateGroup2(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    Trgsw tmp1{param}, tmp2{param}, tmp3{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    auto& pool = ThreadPool::instance();
    for (auto i = 0; i < param.n; i = i + 2) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto& bsk1 = bsk[j];
        auto bsk2 = bsk[j + 1];
        auto bsk3 = bsk[j + 2];
        auto bsk4 = bsk[j + 3];
        auto a12 = a1 + a2;
        auto future1 = pool.enqueue([&]{ rotateTrgsw(bsk2, a2, param); });
        auto future2 = pool.enqueue([&]{ rotateTrgsw(bsk3, a1, param); });
        auto future3 = pool.enqueue([&]{ rotateTrgsw(bsk4, a12, param); });
        future1.get();
        future2.get();
        future3.get();

        for (size_t l = 0; l < param.l; l++) {
            for (size_t k = 0; k < param.k + 1; k++) {
                for (auto k2 = 0; k2 < param.k; k2++) {
                    auto& coeffA1 = bsk1.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA2 = bsk2.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA3 = bsk3.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA4 = bsk4.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffAT = tmp3.trlweSamples[l][k].a[k2].coeffs;
                    for (int n = 0; n < param.N; n++) {
                        coeffAT[n] = coeffA1[n] + coeffA2[n] + coeffA3[n] + coeffA4[n];
                    }
                }
                auto& coeffB1 = bsk1.trlweSamples[l][k].b.coeffs;
                auto& coeffB2 = bsk2.trlweSamples[l][k].b.coeffs;
                auto& coeffB3 = bsk3.trlweSamples[l][k].b.coeffs;
                auto& coeffB4 = bsk4.trlweSamples[l][k].b.coeffs;
                auto& coeffBT = tmp3.trlweSamples[l][k].b.coeffs;
                for (int n = 0; n < param.N; n++) {
                    coeffBT[n] = coeffB1[n] + coeffB2[n] + coeffB3[n] + coeffB4[n];
                }
            }
        }
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
    auto& pool = ThreadPool::instance();
    for (auto i = 0; i < param.n; i = i + 2) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto& bsk1 = bskDft[j];
        auto bsk2 = bskDft[j + 1];
        auto bsk3 = bskDft[j + 2];
        auto bsk4 = bskDft[j + 3];
        auto a12 = a1 + a2;
        auto future1 = pool.enqueue([&]{ rotateTrgswNtt(bsk2, a2, param); });
        auto future2 = pool.enqueue([&]{ rotateTrgswNtt(bsk3, a1, param); });
        auto future3 = pool.enqueue([&]{ rotateTrgswNtt(bsk4, a12, param); });
        future1.get();
        future2.get();
        future3.get();

        auto q = NttHexl::getNttHexl().GetModulus();
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
                    uint64_t val1 = AddUIntMod(coeffB1[n], coeffB2[n], q);
                    uint64_t val2 = AddUIntMod(coeffB3[n], coeffB4[n], q);
                    coeffBT[n] = AddUIntMod(val1, val2, q);
                }
            }
        }
        temp = Trlwe{param.k, param.N};
        externalProductTrgswNtt(temp, tmp3, accum, param.lApprox, param);
        accum = std::move(temp);
        j += batchSize;
    }
}

// todo: general n support
void blindRotateGroup3(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    Trgsw tmp{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    auto& pool = ThreadPool::instance();
    for (auto i = 0; i < param.n; i = i + 3) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto a3 = input.a[i + 2];
        auto& bsk1 = bsk[j];
        auto bsk2 = bsk[j + 1];
        auto bsk3 = bsk[j + 2];
        auto bsk4 = bsk[j + 3];
        auto bsk5 = bsk[j + 4];
        auto bsk6 = bsk[j + 5];
        auto bsk7 = bsk[j + 6];
        auto bsk8 = bsk[j + 7];
        auto a12 = a1 + a2;
        auto a13 = a1 + a3;
        auto a23 = a2 + a3;
        auto a123 = a12 + a3;
        auto future1 = pool.enqueue([&]{ rotateTrgsw(bsk2, a3, param); });
        auto future2 = pool.enqueue([&]{ rotateTrgsw(bsk3, a2, param); });
        auto future3 = pool.enqueue([&]{ rotateTrgsw(bsk4, a23, param); });
        auto future4 = pool.enqueue([&]{ rotateTrgsw(bsk5, a1, param); });
        auto future5 = pool.enqueue([&]{ rotateTrgsw(bsk6, a13, param); });
        auto future6 = pool.enqueue([&]{ rotateTrgsw(bsk7, a12, param); });
        auto future7 = pool.enqueue([&]{ rotateTrgsw(bsk8, a123, param); });

        future1.get();
        future2.get();
        future3.get();
        future4.get();
        future5.get();
        future6.get();
        future7.get();

        for (size_t l = 0; l < param.l; l++) {
            for (size_t k = 0; k < param.k + 1; k++) {
                for (auto k2 = 0; k2 < param.k; k2++) {
                    auto& coeffA1 = bsk1.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA2 = bsk2.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA3 = bsk3.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA4 = bsk4.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA5 = bsk5.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA6 = bsk6.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA7 = bsk7.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffA8 = bsk8.trlweSamples[l][k].a[k2].coeffs;
                    auto& coeffAT = tmp.trlweSamples[l][k].a[k2].coeffs;
                    for (int n = 0; n < param.N; n++) {
                        coeffAT[n] = coeffA1[n] + coeffA2[n] + coeffA3[n] + coeffA4[n] + coeffA5[n] + coeffA6[n] + coeffA7[n] + coeffA8[n];
                    }
                }
                auto& coeffB1 = bsk1.trlweSamples[l][k].b.coeffs;
                auto& coeffB2 = bsk2.trlweSamples[l][k].b.coeffs;
                auto& coeffB3 = bsk3.trlweSamples[l][k].b.coeffs;
                auto& coeffB4 = bsk4.trlweSamples[l][k].b.coeffs;
                auto& coeffB5 = bsk5.trlweSamples[l][k].b.coeffs;
                auto& coeffB6 = bsk6.trlweSamples[l][k].b.coeffs;
                auto& coeffB7 = bsk7.trlweSamples[l][k].b.coeffs;
                auto& coeffB8 = bsk8.trlweSamples[l][k].b.coeffs;
                auto& coeffBT = tmp.trlweSamples[l][k].b.coeffs;
                for (int n = 0; n < param.N; n++) {
                    coeffBT[n] = coeffB1[n] + coeffB2[n] + coeffB3[n] + coeffB4[n] + coeffB5[n] + coeffB6[n] + coeffB7[n] + coeffB8[n];
                }
            }
        }
        temp = Trlwe{param.k, param.N};
        externalProductTrgsw(temp, tmp, accum, param);
        accum = std::move(temp);
        j += batchSize;
    }
}

void blindRotateGroup3Ntt(Trlwe& accum, const vector<TrgswDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlwe temp{param.k, param.N};
    TrgswDft tmp{param};
    int j = 0;
    auto batchSize = 1 << param.group;
    auto& pool = ThreadPool::instance();
    for (auto i = 0; i < param.n; i = i + 3) {
        auto a1 = input.a[i];
        auto a2 = input.a[i + 1];
        auto a3 = input.a[i + 2];
        auto& bsk1 = bskDft[j];
        auto bsk2 = bskDft[j + 1];
        auto bsk3 = bskDft[j + 2];
        auto bsk4 = bskDft[j + 3];
        auto bsk5 = bskDft[j + 4];
        auto bsk6 = bskDft[j + 5];
        auto bsk7 = bskDft[j + 6];
        auto bsk8 = bskDft[j + 7];
        auto a12 = a1 + a2;
        auto a13 = a1 + a3;
        auto a23 = a2 + a3;
        auto a123 = a12 + a3;
        auto future1 = pool.enqueue([&]{ rotateTrgswNtt(bsk2, a3, param); });
        auto future2 = pool.enqueue([&]{ rotateTrgswNtt(bsk3, a2, param); });
        auto future3 = pool.enqueue([&]{ rotateTrgswNtt(bsk4, a23, param); });
        auto future4 = pool.enqueue([&]{ rotateTrgswNtt(bsk5, a1, param); });
        auto future5 = pool.enqueue([&]{ rotateTrgswNtt(bsk6, a13, param); });
        auto future6 = pool.enqueue([&]{ rotateTrgswNtt(bsk7, a12, param); });
        auto future7 = pool.enqueue([&]{ rotateTrgswNtt(bsk8, a123, param); });

        future1.get();
        future2.get();
        future3.get();
        future4.get();
        future5.get();
        future6.get();
        future7.get();

        auto q = NttHexl::getNttHexl().GetModulus();
        auto N = param.N;
        std::vector vecs(6, std::vector<uint64_t>(N));
        for (size_t l = 0; l < param.l; l++) {
            for (size_t k = 0; k < param.k + 1; k++) {
                for (auto k2 = 0; k2 < param.k; k2++) {
                    auto& coeffA1 = bsk1.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA2 = bsk2.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA3 = bsk3.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA4 = bsk4.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA5 = bsk5.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA6 = bsk6.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA7 = bsk7.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffA8 = bsk8.trlweDftSamples[l][k].a[k2].coeffs;
                    auto& coeffAT = tmp.trlweDftSamples[l][k].a[k2].coeffs;

                    EltwiseAddMod(vecs[0].data(), coeffA1.data(), coeffA2.data(), N, q);
                    EltwiseAddMod(vecs[1].data(), coeffA3.data(), coeffA4.data(), N, q);
                    EltwiseAddMod(vecs[2].data(), coeffA5.data(), coeffA6.data(), N, q);
                    EltwiseAddMod(vecs[3].data(), coeffA7.data(), coeffA8.data(), N, q);

                    EltwiseAddMod(vecs[4].data(), vecs[0].data(), vecs[1].data(), N, q);
                    EltwiseAddMod(vecs[5].data(), vecs[2].data(), vecs[3].data(), N, q);

                    EltwiseAddMod(coeffAT.data(), vecs[4].data(), vecs[5].data(), N, q);
                }
                auto& coeffB1 = bsk1.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB2 = bsk2.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB3 = bsk3.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB4 = bsk4.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB5 = bsk5.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB6 = bsk6.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB7 = bsk7.trlweDftSamples[l][k].b.coeffs;
                auto& coeffB8 = bsk8.trlweDftSamples[l][k].b.coeffs;
                auto& coeffBT = tmp.trlweDftSamples[l][k].b.coeffs;

                EltwiseAddMod(vecs[0].data(), coeffB1.data(), coeffB2.data(), N, q);
                EltwiseAddMod(vecs[1].data(), coeffB3.data(), coeffB4.data(), N, q);
                EltwiseAddMod(vecs[2].data(), coeffB5.data(), coeffB6.data(), N, q);
                EltwiseAddMod(vecs[3].data(), coeffB7.data(), coeffB8.data(), N, q);

                EltwiseAddMod(vecs[4].data(), vecs[0].data(), vecs[1].data(), N, q);
                EltwiseAddMod(vecs[5].data(), vecs[2].data(), vecs[3].data(), N, q);

                EltwiseAddMod(coeffBT.data(), vecs[4].data(), vecs[5].data(), N, q);
            }
        }
        temp = Trlwe{param.k, param.N};
        externalProductTrgswNtt(temp, tmp, accum, param.lApprox, param);
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
        case 3:
            blindRotateGroup3(accum, bsk, input, param);
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
        case 3:
            blindRotateGroup3Ntt(accum, bskDft, input, param);
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

void blindRotateWithPreRotNtt(Trlwe& accum, const vector<Trlwe>& trlwe, const vector<vector<TrgswMPDft>>& trgsws,
                              const ScaledTlwe& input, const YatfheParameters& param) {
    const auto n = param.n;
#ifdef TERNARY
    {
        Trlwe tmp{param}, tmp2{param};
        rotateTrlwe(tmp, trlwe[0], input.a[0]);
        rotateTrlwe(tmp2, trlwe[1], -input.a[0]);
        addTrlwe(tmp, tmp, tmp2, trlwe[2]);
        rotateTrlwe(accum, tmp, -input.b);
    }

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        const auto ai = input.a[i];
        const int j = i - 1;
        auto keyOut = trgsws[j][0];
        auto keyOut2 = trgsws[j][1];
        if (ai != 0) {
            rotateTrgswMPNtt(keyOut, ai, param);
            rotateTrgswMPNtt(keyOut2, -ai, param);
        }
        addTrgswMPNtt(keyOut, keyOut, keyOut2, trgsws[j][2]);
        externalProductTrgswMPNttInPlace(accum, keyOut, param.lApprox, param);
    }
#else
    {
        Trlwe tmp{param.k, param.N};
        rotateTrlwe(tmp, trlwe[0], input.a[0]);
        addTrlwe(tmp, tmp, trlwe[1]);
        rotateTrlwe(accum, tmp, -input.b);
    }

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        const auto ai = input.a[i];
        const int j = i - 1;
        auto keyOut = trgsws[j][0];
        if (ai != 0) {
            rotateTrgswMPNtt(keyOut, ai, param);
        }
        addTrgswMPNtt(keyOut, keyOut, trgsws[j][1]);
        externalProductTrgswMPNttInPlace(accum, keyOut, param.lApprox, param);
    }
#endif
}

void blindRotateWithPreRotNttMT(Trlwe& accum, const vector<Trlwe>& trlwe, const vector<vector<TrgswMPDft>>& trgsws,
                                const ScaledTlwe& input, const YatfheParameters& param) {
    const auto n = param.n;
    vector trgswMPDft(param.n-1, TrgswMPDft{param});
#ifdef TERNARY
    preRotateTernary2(accum, trgswMPDft, trlwe, trgsws, input, param.batchSize, param.tasksPerThread, param);

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        externalProductTrgswMPNttInPlace(accum, trgswMPDft[i-1], param.lApprox, param);
    }
#else
    preRotateBinary(accum, trgswMPDft, trlwe, trgsws, input, param.batchSize, param);

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        externalProductTrgswMPNttInPlace(accum, trgswMPDft[i-1], param.lApprox, param);
    }
#endif
}

void blindRotateJP22Ntt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
#ifdef TERNARY
    const auto n = param.n;
    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        auto key0 = bskDft[i][0];
        auto key1 = bskDft[i][1];
        rotateTrgswMPMinusOneNtt(key0, input.a[i], param);
        rotateTrgswMPMinusOneNtt(key1, -input.a[i], param);
        addTrgswMPNtt(key0, key0, key1);

        externalProductTrgswMPNtt(tmp, key0, accum, param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateJP22NttMT(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
#ifdef TERNARY
    const auto n = param.n;
    vector rotated(n, TrgswMPDft{param});

    // pre-rotation in parallel
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(param.batchSize);
    for (int start = 0; start < n; start += param.batchSize) {
        futures.clear();
        const int end = min(start + param.batchSize, n);
        for (int i = start; i < end; ++i) {
            const auto ai = input.a[i];
            futures.emplace_back(pool.enqueue([i, ai, &param, &rotated, &bskDft] {
                Trlwe tmp{param};
                auto key0 = bskDft[i][0];
                auto key1 = bskDft[i][1];
                rotateTrgswMPMinusOneNtt(key0, ai, param);
                rotateTrgswMPMinusOneNtt(key1, -ai, param);
                addTrgswMPNtt(rotated[i], key0, key1);
            }));
        }
        for (auto& f : futures) {
            f.get();
        }
    }

    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        externalProductTrgswMPNtt(tmp, rotated[i], accum, param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateMP21Ntt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
#ifdef TERNARY
    const auto n = param.n;
    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], param.lApprox, param);
        accumulateTrlwe(accum, tmp);

        rotateTrlweMinusOne(tmp, accum, -input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][1], param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], param.lApprox, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateMPInternalNtt(TrgswMP& accum, const vector<TrgswMPDft>& trgsws, const ScaledTlwe& input, const YatfheParameters& param) {
    TrgswMP temp{param};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = accum;
        TrgswMP tmp {param};
        rotateTrgswMP(accum, input.a[i], param);
        subTrgswMP(tmp, accum, temp);
        internalProductTrgswMPNtt(accum, tmp, trgsws[i], param.lApprox, param); // res *= bskI
        addTrgswMP(tmp, accum, temp); // res += input
        accum = std::move(tmp);
    }
}

void blindRotateExternalGeneralNtt(Trlev& accum, const vector<TrgswMPDft>& trgsws, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlev temp{param};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = accum;
        Trlev tmp{param};
        rotateTrlev(accum, input.a[i], param);
        subTrlev(tmp, accum, temp);
        generalExternalProductTrgswMPNtt(accum, trgsws[i], accum, param.lApprox, param);
        addTrlev(tmp, accum, temp); // res += input
        accum = std::move(tmp);
    }
}

void blindRotateExternalGeneralPreRotNtt(Trlev& accum, const vector<TrgswMPDft>& trgsws, const ScaledTlwe& input, const YatfheParameters& param) {
    Trlev temp{param};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        temp = accum;
        Trlev tmp{param};
        rotateTrlev(accum, input.a[i], param);
        subTrlev(tmp, accum, temp);
        generalExternalProductTrgswMPNtt(accum, trgsws[i], accum, param.lApprox, param);
        addTrlev(tmp, accum, temp); // res += input
        accum = std::move(tmp);
    }
}

void blindRotateInternalPairWiseNtt(Trlwe& accum, vector<vector<TrgswMP>>& trgsws, vector<vector<TrgswMPDft>>& trgswDfts,
                                    const ScaledTlwe& input, const int batchSize, const YatfheParameters& param) {
    auto n = param.n;
    vector trgswMP(n / 2, TrgswMP{param});
    vector trgswMPDft(n / 2, TrgswMPDft{param});
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    preRotateInternal(trgswMP, trgswMPDft, trgsws, trgswDfts, input.a, param);
    while (n > 1) {
        const auto newSize = n / 2;
        vector newTrgswMP(newSize, TrgswMP{param});
        vector newTrgswMPDft(newSize, TrgswMPDft{param});

//         for (int i = 0; i < n / 2; i++) {
//             if (i % 2 == 0) {
//                 internalProductTrgswMPNtt(newTrgswMP[i/2], trgswMP[i], trgswMPDft[i], param);
//             } else {
//                 internalProductTrgswMPNtt(newTrgswMPDft[(i-1)/2], trgswMP[i], trgswMPDft[i], param);
//             }
//         }

        for (int start = 0; start < n / 2; start += batchSize) {
            futures.clear();
            for (int i = start; i < min(start + batchSize, n / 2); i++) {
                futures.emplace_back(pool.enqueue([i, &newTrgswMP, &newTrgswMPDft, &trgswMP, &trgswMPDft, &param] {
                   if (i % 2 == 0) {
                       internalProductTrgswMPNtt(newTrgswMPDft[i/2], trgswMP[i], trgswMPDft[i], param.l, param);
                   } else {
                       internalProductTrgswMPNtt(newTrgswMP[(i-1)/2], trgswMP[i], trgswMPDft[i], param.l, param);
                   }
                }));
            }
            for (auto& f : futures) {
                f.wait();
            }
        }
        trgswMP = std::move(newTrgswMP);
        trgswMPDft = std::move(newTrgswMPDft);
        n = newSize;
    }
    Trlwe cop = accum;
    externalProductTrgswMPNtt(accum, trgswMPDft[0], cop, param.lApprox, param);
}

void blindRotateInternalPairWiseAsymNtt(Trlwe& accum, vector<vector<Trlev>>& trlevs, vector<vector<TrgswMPDft>>& trgswDfts,
    const ScaledTlwe& input, const TrlevDft& sSquare, const int batchSize, const YatfheParameters& param) {
    auto n = param.n;
    vector trlev(n / 2, Trlev{param});
    vector trgswMPDft(n / 2, TrgswMPDft{param});
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    preRotateInternalAsym(trlev, trgswMPDft, trlevs, trgswDfts, input.a, param);
    while (n > 1) {
        const auto newSize = n / 2;
        vector newTrlev(newSize, Trlev{param});
        vector newTrgswMPDft(newSize, TrgswMPDft{param});
        for (int start = 0; start < newSize; start += batchSize) {
            futures.clear();
            for (int i = start; i < min(start + batchSize, newSize); i++) {
                futures.emplace_back(pool.enqueue([i, &newTrlev, &newTrgswMPDft, &trlev, &trgswMPDft, &sSquare, &param] {
                   if (i % 2 == 0) {
                       internalProductAsymTrgswMPNtt(newTrgswMPDft[i/2], trgswMPDft[i], trlev[i], sSquare, param.l, param);
                   } else {
                       generalExternalProductTrgswMPNtt(newTrlev[(i-1)/2], trgswMPDft[i], trlev[i], param.l, param);
                   }
                }));
            }
            for (auto& f : futures) {
                f.wait();
            }
        }
        trlev = std::move(newTrlev);
        trgswMPDft = std::move(newTrgswMPDft);
        n = newSize;
    }
    Trlwe copy = accum;
    externalProductTrgswMPNtt(accum, trgswMPDft[0], copy, param.lApprox, param);
}

void blindRotateInternalPairWiseAsymOptNtt(Trlwe& out, vector<vector<Trlev>>& trlevs, vector<Trlwe>& trlwes,
                                           vector<vector<TrgswMPDft>>& trgswDfts, const ScaledTlwe& input,
                                           const TrlevDft& sSquare, const int batchSize, const YatfheParameters& param) {
    auto n = param.n;
    vector trlev(n / 2, Trlev{param});
    vector trgswMPDft(n / 2, TrgswMPDft{param});
    Trlwe trlwe{param.k, param.N};
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    preRotateInternalAsymOpt(trlev, trlwe, trgswMPDft, trlevs, trlwes, trgswDfts, input.a, param);
    while (n > 1) {
        const auto newSize = n / 2;
        vector newTrlev(newSize, Trlev{param});
        vector newTrgswMPDft(newSize, TrgswMPDft{param});
        Trlwe newTrlwe {param.k, param.N};
        for (int start = 0; start < newSize; start += batchSize) {
            futures.clear();
            for (int i = start; i < min(start + batchSize, newSize); i++) {
                futures.emplace_back(pool.enqueue([i, newSize, &newTrlev, &newTrlwe, &newTrgswMPDft, &trlev, &trlwe, &trgswMPDft, &sSquare, &param] {
                    if (i == newSize - 1) {
                        externalProductTrgswMPNtt(newTrlwe, trgswMPDft[i], trlwe, param.lApprox, param);
                        return;
                    }
                    if (i % 2 == 0) {
                        internalProductAsymTrgswMPNtt(newTrgswMPDft[i/2], trgswMPDft[i], trlev[i], sSquare, param.l, param);
                    } else {
                        generalExternalProductTrgswMPNtt(newTrlev[(i-1)/2], trgswMPDft[i], trlev[i], param.l, param);
                    }
                }));
            }
            for (auto& f : futures) {
                f.wait();
            }
        }
        trlev = std::move(newTrlev);
        trlwe = std::move(newTrlwe);
        trgswMPDft = std::move(newTrgswMPDft);
        n = newSize;
    }
    rotateTrlwe(out, trlwe, -input.b);
}
