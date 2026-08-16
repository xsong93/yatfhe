//
// Created by xintong on 4/21/25.
//
#include <stdexcept>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/cmux.h"
#include "yatfhe/ntt_hexl.h"
#include "yautil/time_counter.h"
#include "yautil/multi_threading.h"
#include "yautil/ya_serializer.h"

namespace {
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
                    if (ai != 0) {
                        rotateTrgswMPNtt(keyOut, keyIn[0], ai, param);
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

                        TrgswMPDft keyOut2{param};
                        if (ai != 0) {
                            rotateTrgswMPNtt(keyOut, keyIn[0], ai, param);
                            rotateTrgswMPNtt(keyOut2, keyIn[1], -ai, param);
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

    void switchSchemeInBatchBinary(vector<vector<TrgswMPDft>>& bsk, const TrlevDft& s2, const ScaledTlwe& input,
                                   const YatfheParameters& param) {
        /**
         * Performs:
         *  for (auto i = 0; i < n-1; i++) {
         *       bsk[i][0].c = vector(level, vector(param.k, TrlweDft(param.k, param.N)));
         *       for (auto l = 0; l < level; l++) {
         *           switchTrlweToSecretEmbeddingNtt(bsk[i][0].c[l], bsk[i][0].cPrime[l], s2, param);
         *       }
         *   }
         */
        const auto n = param.n;
        const auto batchSize = param.batchSize;
        const auto tasksPerThread = param.tasksPerThread;
        const auto level = bsk[0][0].l;
        auto& pool = ThreadPool::instance();

        vector<future<void>> futures;
        futures.reserve(batchSize);

        for (int start = 0; start < n - 1; start += batchSize * tasksPerThread) {
            futures.clear();
            const int end = min(start + batchSize * tasksPerThread, n - 1);

            // Process tasks in chunks of 'tasksPerThread'
            for (int chunkStart = start; chunkStart < end; chunkStart += tasksPerThread) {
                const int chunkEnd = min(chunkStart + tasksPerThread, end);

                futures.emplace_back(pool.enqueue([&bsk, &s2, &input, &param, chunkStart, chunkEnd, level] {
                    for (int i = chunkStart; i < chunkEnd; ++i) {
                        if (input.a[i+1] == 0) {
                            continue;
                        }
                        bsk[i][0].c.resize(level, vector(param.k, TrlweDft(param.k, param.N)));
                        for (auto l = 0; l < level; l++) {
                            switchTrlweToSecretEmbeddingNtt(bsk[i][0].c[l], bsk[i][0].cPrime[l], s2, param);
                        }
                    }
                }));
            }
            for (auto& f : futures) {
                f.get();
            }
        }
    }

    void switchSchemeInBatchBinaryOpt(vector<vector<TrgswMPDft>>& bsk, const TrlevDft& s2, const ScaledTlwe& input,
                                   const YatfheParameters& param,
                                   const vector<vector<vector<DecompPolynomial>>>& bskDecompA) {
        const auto n = param.n;
        const auto batchSize = param.batchSize;
        const auto tasksPerThread = param.tasksPerThread;
        const auto level = bsk[0][0].l;
        auto& pool = ThreadPool::instance();

        vector<future<void>> futures;
        futures.reserve(batchSize);

        for (int start = 0; start < n - 1; start += batchSize * tasksPerThread) {
            futures.clear();
            const int end = min(start + batchSize * tasksPerThread, n - 1);

            // Process tasks in chunks of 'tasksPerThread'
            for (int chunkStart = start; chunkStart < end; chunkStart += tasksPerThread) {
                const int chunkEnd = min(chunkStart + tasksPerThread, end);

                futures.emplace_back(pool.enqueue([&bsk, &s2, &input, &param, &bskDecompA, chunkStart, chunkEnd, level] {
                    for (int i = chunkStart; i < chunkEnd; ++i) {
                        if (input.a[i+1] == 0) {
                            continue;
                        }
                        bsk[i][0].c.resize(level);
                        for (auto l = 0; l < level; l++) {
                            auto& c = bsk[i][0].c[l];
                            auto& cPrime = bsk[i][0].cPrime[l];
                            auto& decompA = bskDecompA[i * level + l];
                            c.resize(param.k, TrlweDft(param.k, param.N));
                            cPrime.a.resize(param.k, NttPolynomial(param.N));
                            switchTrlweToSecretEmbeddingNttOpt(bsk[i][0].c[l], bsk[i][0].cPrime[l],
                                decompA, s2, param);
                        }
                    }
                }));
            }
            for (auto& f : futures) {
                f.get();
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

    void rotateBskComponent(const vector<vector<DecompPolynomial>>& bskDecompB, const vector<vector<vector<DecompPolynomial>>>& bskDecompA,
                            vector<vector<vector<DecompPolynomial>>>& rotatedA, vector<vector<DecompPolynomial>>& rotatedB,
                            const int32_t a, const int keyIndex, const int level, const YatfheParameters& param) {
        const auto q = NttHexl::getNttHexl().GetModulus();
        for (auto lvl = 0; lvl < level; lvl++) {
            const auto& decompA = bskDecompA[keyIndex * level + lvl];
            const auto& decompB = bskDecompB[keyIndex * level + lvl];
            for (auto dl = 0; dl < param.l; dl++) {
                for (auto k1 = 0; k1 < param.k; k1++) {
                    rotateDecompPolynomialMinusOne(rotatedA[lvl][dl][k1], a, decompA[dl][k1]);
                }
                rotateDecompPolynomialMinusOne(rotatedB[lvl][dl], a, decompB[dl]);

                // b+1
                if (lvl == dl) rotatedB[lvl][lvl].coeffs[0] = static_cast<Decomp>(rotatedB[lvl][lvl].coeffs[0] + 1);
            }
        }
    }

    void deserializeAndRotateBskComponent(vector<vector<DecompPolynomial>>& bskDecompB, vector<vector<vector<DecompPolynomial>>>& bskDecompA,
                                       vector<vector<vector<DecompPolynomial>>>& rotatedA, vector<vector<DecompPolynomial>>& rotatedB,
                                       const int32_t a, const int keyIndex, const int level, std::ifstream& inFile,
                                       const YatfheParameters& param) {
        for (auto l = 0; l < level; l++) {
            deserializeNestedVector(bskDecompA[keyIndex * level + l], inFile);
            deserializeNestedVector(bskDecompB[keyIndex * level + l], inFile);
        }

        for (auto lvl = 0; lvl < level; lvl++) {
            const auto& decompA = bskDecompA[keyIndex * level + lvl];
            const auto& decompB = bskDecompB[keyIndex * level + lvl];
            for (auto dl = 0; dl < param.l; dl++) {
                for (auto k1 = 0; k1 < param.k; k1++) {
                    rotateDecompPolynomialMinusOne(rotatedA[lvl][dl][k1], a, decompA[dl][k1]);
                }
                rotateDecompPolynomialMinusOne(rotatedB[lvl][dl], a, decompB[dl]);

                // b+1
                if (lvl == dl) rotatedB[lvl][lvl].coeffs[0] = static_cast<Decomp>(rotatedB[lvl][lvl].coeffs[0] + 1);
            }
        }
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
        TrgswMPDft keyOut{param};
        TrgswMPDft keyOut2{param};
        if (ai != 0) {
            rotateTrgswMPNtt(keyOut, trgsws[j][0], ai, param);
            rotateTrgswMPNtt(keyOut2, trgsws[j][1], -ai, param);
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
        TrgswMPDft keyOut{param};
        if (ai != 0) {
            rotateTrgswMPNtt(keyOut, trgsws[j][0], ai, param);
        }
        addTrgswMPNtt(keyOut, keyOut, trgsws[j][1]);
        externalProductTrgswMPNttInPlace(accum, keyOut, param.lApprox, param);
    }
#endif
}

void blindRotateWithPreRotNttMT(Trlwe& accum, const vector<Trlwe>& trlwe, const vector<vector<TrgswMPDft>>& trgsws,
                                const ScaledTlwe& input, const YatfheParameters& param) {
    const auto n = param.n;
    const auto level = trgsws[0][0].l;
    vector trgswMPDft(param.n-1, TrgswMPDft{param, level});
#ifdef TERNARY
    preRotateTernary(accum, trgswMPDft, trlwe, trgsws, input, param.batchSize, param.tasksPerThread, param);

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        externalProductTrgswMPNttInPlace(accum, trgswMPDft[i-1], level, param);
    }
#else
    preRotateBinary(accum, trgswMPDft, trlwe, trgsws, input, param.batchSize, param);

    for (auto i = 1; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        externalProductTrgswMPNttInPlace(accum, trgswMPDft[i-1], level, param);
    }
#endif
}

void blindRotateOptNtt(Trlwe& accum, const vector<Trlwe>& bskFirst,  const vector<vector<TrgswMPDft>>& bskDft,
                       const ScaledTlwe& input, const TorusPolynomial& v, const YatfheParameters& param) {
    const auto level = bskDft[0][0].l;
    const auto n = param.n;

#ifdef TERNARY
    // handle first RLwe key component
    // R(v) + (X^a0 - 1)R(v*s00) + (X^{-a0} - 1)R(v*s01)
    {
        Trlwe tmp{param}, tmp2{param};
        rotateTrlweMinusOne(tmp, bskFirst[0], input.a[0]);
        rotateTrlweMinusOne(tmp2, bskFirst[1], -input.a[0]);
        addTrlwe(tmp, tmp, tmp2);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);
    }

    // accumulate on the remaining n-1 key components
    for (auto i = 0; i < n-1; i++) {
        const int j = i + 1;
        if (input.a[j] == 0) {
            continue;
        }
        Trlwe tmp{param};
        TrgswMPDft key0{param};
        TrgswMPDft key1{param};
        rotateTrgswMPMinusOneNtt(key0, bskDft[i][0], input.a[j], param);
        rotateTrgswMPMinusOneNtt(key1, bskDft[i][1], -input.a[j], param);
        addTrgswMPNtt(key0, key0, key1);

        externalProductTrgswMPNtt(tmp, key0, accum, level, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    // handle first Rlwe key component
    // R(v) + (X^a0 - 1)R(v*s0)
    {
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);
    }

    // accumulate on the remaining n-1 key components
    for (auto i = 0; i < n-1; i++) {
        if (input.a[i+1] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i+1]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

//lazy init
void blindRotateLazyNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk, bool& initialized,
                        const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2, const YatfheParameters& param) {
    const auto level = bsk[0][0].l;
    const auto n = param.n;

#ifdef TERNARY
    // handle first Rlwe key component
    // R(v) + (X^a0 - 1)R(v*s00) + (X^{-a0} - 1)R(v*s01)
    {
        Trlwe tmp{param}, tmp2{param};
        rotateTrlweMinusOne(tmp, bskFirst[0], input.a[0]);
        rotateTrlweMinusOne(tmp2, bskFirst[1], -input.a[0]);
        addTrlwe(tmp, tmp, tmp2);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);
    }

    // accumulate on the remaining n-1 key components
    for (auto i = 0; i < n-1; i++) {
        const int j = i + 1;
        if (input.a[j] == 0) {
            continue;
        }
        Trlwe tmp{param};
        TrgswMPDft key0{param};
        TrgswMPDft key1{param};
        rotateTrgswMPMinusOneNtt(key0, bskDft[i][0], input.a[j], param);
        rotateTrgswMPMinusOneNtt(key1, bskDft[i][1], -input.a[j], param);
        addTrgswMPNtt(key0, key0, key1);

        externalProductTrgswMPNtt(tmp, key0, accum, level, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    // handle first Rlwe key component
    // R(v) + (X^a0 - 1)R(v*s0)
    {
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);
    }

    // calculate secret dependent part of the bsk
    if (!initialized) {
        switchSchemeInBatchBinary(bsk, s2, input, param);
        initialized = true;
    }

    // accumulate on the remaining n-1 key components
    for (auto i = 0; i < n-1; i++) {
        if (input.a[i+1] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i+1]);
        externalProductTrgswMPNttInPlace(tmp, bsk[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateLazyMTNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk,
                           const vector<vector<vector<DecompPolynomial>>>& bskDecompA,
                           const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2,
                           const YatfheParameters& param) {
    const auto level = bsk[0][0].l;
    const auto n = param.n;

#ifdef TERNARY
    throw std::runtime_error("blindRotateLazyMTNtt: not implemented for ternary keys");
#else
    // handle first Rlwe key component
    // R(v) + (X^a0 - 1)R(v*s0)
    {
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);
    }

    // calculate secret dependent part of the bsk
    switchSchemeInBatchBinaryOpt(bsk, s2, input, param, bskDecompA);

    // accumulate on the remaining n-1 key components
    for (auto i = 0; i < n-1; i++) {
        if (input.a[i+1] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i+1]);
        externalProductTrgswMPNttInPlace(tmp, bsk[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateLazyPipeNtt(Trlwe& accum, const BootstrappingKeyMPLazyPipe& bsk, const ScaledTlwe& input,
                            const TorusPolynomial& v, const YatfheParameters& param) {
    const auto level = param.lApprox;
    const auto n = param.n;
    TrgswMPDft rotated0{param, level};
    TrgswMPDft rotated1{param, level};
    vector rotatedADecomp0(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector rotatedADecomp1(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector rotatedBDecomp0(level, vector(param.l, DecompPolynomial{param.N}));
    vector rotatedBDecomp1(level, vector(param.l, DecompPolynomial{param.N}));
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(1 + level);

#ifdef TERNARY
    throw std::runtime_error("blindRotateLazyPipeNtt: not implemented for ternary keys");
#else
    auto& s2 = bsk.s2Dft;
    auto& bskDecompA = bsk.bskDecompA;
    auto& bskDecompB = bsk.bskDecompB;

    // handle first two key components
    {
        const auto a1 = input.a[1];
        futures.emplace_back(pool.enqueue([&bskDecompB, &bskDecompA, &rotatedADecomp0, &rotatedBDecomp0, a1, level, &param] {
            rotateBskComponent(bskDecompB, bskDecompA, rotatedADecomp0, rotatedBDecomp0, a1, 0, level, param);
        }));

        // prep
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bsk.bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);

        for (auto& f : futures) f.get();
        futures.clear();
    }

    // main pipeline
    for (auto i = 0; i < n; i++) {
        auto& currRotated = i % 2 == 0 ? rotated0 : rotated1;
        auto& nextRotated = i % 2 == 0 ? rotated1 : rotated0;
        auto& currRotatedA = i % 2 == 0 ? rotatedADecomp0 : rotatedADecomp1;
        auto& nextRotatedA = i % 2 == 0 ? rotatedADecomp1 : rotatedADecomp0;
        auto& currRotatedB = i % 2 == 0 ? rotatedBDecomp0 : rotatedBDecomp1;
        auto& nextRotatedB = i % 2 == 0 ? rotatedBDecomp1 : rotatedBDecomp0;

        // automorphism
        if (i < n - 2) {
            const auto aNext = input.a[i + 2];
            futures.emplace_back(pool.enqueue([&bskDecompB, &bskDecompA, &nextRotatedA, &nextRotatedB, aNext, level, i, &param] {
                rotateBskComponent(bskDecompB, bskDecompA, nextRotatedA, nextRotatedB, aNext, i + 1, level, param);
            }));
        }

        // scheme switching
        if (i < n - 1) {
            for (auto l = 0; l < level; l++) {
                futures.emplace_back(pool.enqueue([&nextRotated, &currRotatedA, &currRotatedB, &s2, l, &param] {
                    clearTrlwe(nextRotated.cPrime[l]);
                    for (auto& item : nextRotated.c[l]) {
                        clearTrlwe(item);
                    }
                    switchDecompTrlweToSecretEmbeddingNtt(nextRotated.c[l], nextRotated.cPrime[l], currRotatedA[l], currRotatedB[l], s2, param);
                }));
            }
        }

        // main thread: accumulation
        if (i >= 1 && input.a[i] != 0) {
            externalProductTrgswMPNttInPlace(accum, currRotated, level, param);
        }
        if (i < n - 1) {
            for (auto& f : futures) f.get();
            futures.clear();
        }
    }
#endif
}



void blindRotatePipeInitNtt(Trlwe& accum, BootstrappingKeyMPLazyPipe& bsk, const ScaledTlwe& input, const TorusPolynomial& v,
                                   const string& fileName, const YatfheParameters& param) {
    int level;
    std::ifstream inFile(fileName, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");
    readPOD(inFile, level);
    const auto n = param.n;
    TrgswMPDft rotated0{param, level};
    TrgswMPDft rotated1{param, level};
    vector rotatedADecomp0(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector rotatedADecomp1(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector rotatedBDecomp0(level, vector(param.l, DecompPolynomial{param.N}));
    vector rotatedBDecomp1(level, vector(param.l, DecompPolynomial{param.N}));
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(1 + level);


#ifdef TERNARY
    throw std::runtime_error("blindRotatePipeInitNtt: not implemented for ternary keys");
#else

    // read keys
    bsk.bskFirst.resize(1);
    bsk.bskDecompA.resize((n - 1) * level);
    bsk.bskDecompB.resize((n - 1) * level);
    deserialize(bsk.bskFirst[0], inFile);
    deserialize(bsk.s2Dft, inFile);
    auto& s2 = bsk.s2Dft;
    auto& bskDecompA = bsk.bskDecompA;
    auto& bskDecompB = bsk.bskDecompB;

    // handle first two key components
    {
        const auto a1 = input.a[1];
        futures.emplace_back(pool.enqueue([&bskDecompB, &bskDecompA, &rotatedADecomp0, &rotatedBDecomp0, a1, level, &inFile, &param] {
            deserializeAndRotateBskComponent(bskDecompB, bskDecompA, rotatedADecomp0, rotatedBDecomp0, a1, 0, level, inFile, param);
        }));

        // prep
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bsk.bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);

        for (auto& f : futures) f.get();
        futures.clear();
    }

    // main pipeline
    for (auto i = 0; i < n; i++) {
        auto& currRotated = i % 2 == 0 ? rotated0 : rotated1;
        auto& nextRotated = i % 2 == 0 ? rotated1 : rotated0;
        auto& currRotatedA = i % 2 == 0 ? rotatedADecomp0 : rotatedADecomp1;
        auto& nextRotatedA = i % 2 == 0 ? rotatedADecomp1 : rotatedADecomp0;
        auto& currRotatedB = i % 2 == 0 ? rotatedBDecomp0 : rotatedBDecomp1;
        auto& nextRotatedB = i % 2 == 0 ? rotatedBDecomp1 : rotatedBDecomp0;

        // automorphism
        if (i < n - 2) {
            const auto aNext = input.a[i + 2];
            futures.emplace_back(pool.enqueue([&bskDecompB, &bskDecompA, &nextRotatedA, &nextRotatedB, aNext, level, i, &inFile, &param] {
                deserializeAndRotateBskComponent(bskDecompB, bskDecompA, nextRotatedA, nextRotatedB, aNext, i + 1, level, inFile, param);
            }));
        }

        // scheme switching
        if (i < n - 1) {
            for (auto l = 0; l < level; l++) {
                futures.emplace_back(pool.enqueue([&nextRotated, &currRotatedA, &currRotatedB, &s2, l, &param] {
                    clearTrlwe(nextRotated.cPrime[l]);
                    for (auto& item : nextRotated.c[l]) {
                        clearTrlwe(item);
                    }
                    switchDecompTrlweToSecretEmbeddingNtt(nextRotated.c[l], nextRotated.cPrime[l], currRotatedA[l], currRotatedB[l], s2, param);
                }));
            }
        }

        // main thread: accumulation
        if (i >= 1 && input.a[i] != 0) {
            externalProductTrgswMPNttInPlace(accum, currRotated, level, param);
        }
        if (i < n - 1) {
            for (auto& f : futures) f.get();
            futures.clear();
        }
    }
    inFile.close();
#endif
}

void blindRotateLazyPipeAltNtt(Trlwe& accum, const BootstrappingKeyMPLazyPipeAlt& bsk, const ScaledTlwe& input, const TorusPolynomial& v,
                               const YatfheParameters& param) {
    const auto level = param.lApprox;
    const auto n = param.n;
    TrgswMPDft expanded0{param, level};
    TrgswMPDft expanded1{param, level};
    vector decompA0(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector decompA1(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector b0(level, TorusPolynomial{param.N});
    vector b1(level, TorusPolynomial{param.N});
    auto& pool = ThreadPool::instance();
    TaskGroup group;

#ifdef TERNARY
    throw std::runtime_error("blindRotateLazyPipeAltNtt: not implemented for ternary keys");
#else
    auto& s2 = bsk.s2Dft;
    auto& bskPrime = bsk.bskPrime;
    vector<Trlwe> holders(level, Trlwe{param});

    int keyIndex = 0;
    int rotateBy = input.a[1];
    auto* currExpanded = &expanded0;
    auto* nextExpanded = &expanded1;
    auto* currDecompA = &decompA0;
    auto* nextDecompA = &decompA0;
    auto* currB = &b0;
    auto* nextB = &b0;

    // automorphism
    const auto stageA = [&](const int l) {
        rotateTrlweMinusOneBPlusOne(holders[l], (*nextB)[l], bskPrime[keyIndex].cPrime[l], rotateBy,
                                static_cast<Torus>(1) << (param.torusBits - (l + 1) * param.radixBits));
        gadgetDecomposeTrlweA((*nextDecompA)[l], holders[l].a, param);
    };

    // scheme switching
    const auto stageB = [&](const int l) {
        clearTrlwe(nextExpanded->cPrime[l]);
        for (auto &item: nextExpanded->c[l]) {
            clearTrlwe(item);
        }
        switchTrlweToSecretEmbeddingNttMix(nextExpanded->c[l], nextExpanded->cPrime[l], (*currDecompA)[l],
                                           (*currB)[l], s2, param);
    };

    // handle first two key components
    {
        pool.run(group, level, stageA);   // keyIndex 0 into decompA0/b0

        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bsk.bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);

        group.wait();
    }

    // accumulate on the n - 1 key components
    for (auto i = 0; i < n; i++) {
        const bool even = i % 2 == 0;
        currExpanded = even ? &expanded0 : &expanded1;
        nextExpanded = even ? &expanded1 : &expanded0;
        currDecompA = even ? &decompA0 : &decompA1;
        nextDecompA = even ? &decompA1 : &decompA0;
        currB = even ? &b0 : &b1;
        nextB = even ? &b1 : &b0;

        if (i < n - 1) {
            pool.run(group, level, stageB);
        }

        if (i < n - 2) {
            keyIndex = i + 1;
            rotateBy = input.a[i + 2];
            pool.run(group, level, stageA);
        }

        // accumulation
        if (i >= 1 && input.a[i] != 0) {
            externalProductTrgswMPNttInPlace(accum, *currExpanded, level, param);
        }

        if (i < n - 1) {
            group.wait();
        }
    }
#endif
}

void blindRotateLazyPipeAltInitNtt(Trlwe& accum, BootstrappingKeyMPLazyPipeAlt& bsk, const ScaledTlwe& input, const TorusPolynomial& v,
                                   const string& fileName, const YatfheParameters& param) {
    const auto level = param.lApprox;
    const auto n = param.n;
    TrgswMPDft expanded0{param, level};
    TrgswMPDft expanded1{param, level};
    vector decompA0(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector decompA1(level, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    vector b0(level, TorusPolynomial{param.N});
    vector b1(level, TorusPolynomial{param.N});
    auto& pool = ThreadPool::instance();
    TaskGroup group;
    std::ifstream inFile(fileName, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

#ifdef TERNARY
    throw std::runtime_error("blindRotateLazyPipeAltInitNtt: not implemented for ternary keys");
#else
    // read keys
    bsk.bskFirst.resize(1);
    bsk.bskPrime.resize(n - 1);
    deserialize(bsk.bskFirst[0], inFile);
    deserialize(bsk.s2Dft, inFile);
    auto& s2 = bsk.s2Dft;
    auto& bskPrime = bsk.bskPrime;
    vector<Trlwe> holders(level, Trlwe{param});

    int keyIndex = 0;
    int readIndex = 0;
    bool readAhead = false;
    int rotateBy = input.a[1];
    auto* currExpanded = &expanded0;
    auto* nextExpanded = &expanded1;
    auto* currDecompA = &decompA0;
    auto* nextDecompA = &decompA0;
    auto* currB = &b0;
    auto* nextB = &b0;

    // automorphism, with the read folded into slot 0
    const auto stageA = [&](const int l) {
        rotateTrlweMinusOneBPlusOne(holders[l], (*nextB)[l], bskPrime[keyIndex].cPrime[l], rotateBy,
                                static_cast<Torus>(1) << (param.torusBits - (l + 1) * param.radixBits));
        gadgetDecomposeTrlweA((*nextDecompA)[l], holders[l].a, param);
        if (l == 0 && readAhead) {
            deserialize(bskPrime[readIndex], inFile);
        }
    };

    // scheme switching
    const auto stageB = [&](const int l) {
        clearTrlwe(nextExpanded->cPrime[l]);
        for (auto &item: nextExpanded->c[l]) {
            clearTrlwe(item);
        }
        switchTrlweToSecretEmbeddingNttMix(nextExpanded->c[l], nextExpanded->cPrime[l], (*currDecompA)[l],
                                           (*currB)[l], s2, param);
    };

    // handle first two key components
    {
        deserialize(bskPrime[0], inFile);
        readIndex = 1;
        readAhead = true;
        pool.run(group, level, stageA);

        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, bsk.bskFirst[0], input.a[0]);
        addTorusPolynomial(tmp.b, tmp.b, v);
        rotateTrlwe(accum, tmp, -input.b);

        group.wait();
    }

    // accumulate on the n - 1 key components
    for (auto i = 0; i < n; i++) {
        const bool even = i % 2 == 0;
        currExpanded = even ? &expanded0 : &expanded1;
        nextExpanded = even ? &expanded1 : &expanded0;
        currDecompA = even ? &decompA0 : &decompA1;
        nextDecompA = even ? &decompA1 : &decompA0;
        currB = even ? &b0 : &b1;
        nextB = even ? &b1 : &b0;

        // Scheme switching
        if (i < n - 1) {
            pool.run(group, level, stageB);
        }

        // automorphism
        if (i < n - 2) {
            keyIndex = i + 1;
            rotateBy = input.a[i + 2];
            readIndex = i + 2;
            readAhead = i < n - 3;
            pool.run(group, level, stageA);
        }

        // accumulation
        if (i >= 1 && input.a[i] != 0) {
            externalProductTrgswMPNttInPlace(accum, *currExpanded, level, param);
        }

        if (i < n - 1) {
            group.wait();
        }
    }
    inFile.close();
#endif
}

// SFBS blind rotation (WWL+24, Algorithm 3)
void blindRotateWWL24Ntt(Trlwe& accum, const BootstrappingKeyWWL24& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
#ifdef TERNARY
    throw std::runtime_error("blindRotateWWL24Ntt: not implemented for ternary keys");
#else
    const auto& bskDft = bsk.bskDft;
    const auto& s2 = bsk.s2Dft;
    const auto level = bskDft[0].l;
    auto& pool = ThreadPool::instance();
    TaskGroup group;
    vector cRows(level, vector(param.k, TrlweDft(param.k, param.N)));

    int keyIndex = 0;
    const auto schemeSwitch = [&](const int l) {
        switchTrlweToSecretEmbeddingNtt(cRows[l], bskDft[keyIndex].cPrime[l], s2, param);
    };

    Trlwe tmp{param};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        keyIndex = i;
        pool.run(group, level, schemeSwitch);
        group.wait();
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductSplitNttInPlace(tmp, cRows, bskDft[i].cPrime, level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

// SFBS blind rotation that performs NTTs on server
void blindRotateWWL24AltNtt(Trlwe& accum, const BootstrappingKeyWWL24Alt& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
#ifdef TERNARY
    throw std::runtime_error("blindRotateWWL24AltNtt: not implemented for ternary keys");
#else
    const auto& trgsws = bsk.trgsws;
    const auto& s2 = bsk.s2Dft;
    const auto level = trgsws[0].l;
    auto& pool = ThreadPool::instance();
    TaskGroup group;
    vector cRows(level, vector(param.k, TrlweDft(param.k, param.N)));
    vector cPrimeRows(level, TrlweDft(param.k, param.N));

    int keyIndex = 0;
    const auto schemeSwitch = [&](const int l) {
        switchTrlweToSecretEmbeddingAltNtt(cRows[l], cPrimeRows[l], trgsws[keyIndex].cPrime[l], s2, param);
    };

    Trlwe tmp{param};
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        keyIndex = i;
        pool.run(group, level, schemeSwitch);
        group.wait();
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductSplitNttInPlace(tmp, cRows, cPrimeRows, level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}


void blindRotateJP22Ntt(Trlwe& accum, const BootstrappingKeyMP& bsk, const ScaledTlwe& input, const YatfheParameters& param) {
    auto& bskDft = bsk.bskDft;
    const auto level = bskDft[0][0].l;
#ifdef TERNARY
    const auto n = param.n;
    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        TrgswMPDft key0{param};
        TrgswMPDft key1{param};
        rotateTrgswMPMinusOneNtt(key0, bskDft[i][0], input.a[i], param);
        rotateTrgswMPMinusOneNtt(key1, bskDft[i][1], -input.a[i], param);
        addTrgswMPNtt(key0, key0, key1);

        externalProductTrgswMPNtt(tmp, key0, accum, level, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateJP22NttMT(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    const auto level = bskDft[0][0].l;
#ifdef TERNARY
    const auto n = param.n;
    const auto tasksPerThread = param.tasksPerThread;
    vector rotated(n, TrgswMPDft{param, level});

    // pre-rotation in parallel
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(param.batchSize);

    for (int start = 0; start < n; start += param.batchSize * tasksPerThread) {
        futures.clear();
        const int end = min(start + param.batchSize * tasksPerThread, n);

        for (int chunkStart = start; chunkStart < end; chunkStart += tasksPerThread) {
            const int chunkEnd = min(chunkStart + tasksPerThread, end);

            futures.emplace_back(pool.enqueue([chunkStart, chunkEnd, &input, &param, &rotated, &bskDft] {
                for (int i = chunkStart; i < chunkEnd; ++i) {
                    const auto ai = input.a[i];
                    Trlwe tmp{param};
                    TrgswMPDft key0{param};
                    TrgswMPDft key1{param};
                    rotateTrgswMPMinusOneNtt(key0, bskDft[i][0], ai, param);
                    rotateTrgswMPMinusOneNtt(key1, bskDft[i][1], -ai, param);
                    addTrgswMPNtt(rotated[i], key0, key1);
                }
            }));
        }
        for (auto& f : futures) {
            f.get();
        }
    }

    // accumulation
    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        externalProductTrgswMPNtt(tmp, rotated[i], accum, level, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        Trlwe tmp{param};
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
}

void blindRotateMP21Ntt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param) {
    const auto level = bskDft[0][0].l;
    Trlwe tmp{param};
#ifdef TERNARY
    const auto n = param.n;
    for (auto i = 0; i < n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], level, param);
        accumulateTrlwe(accum, tmp);

        rotateTrlweMinusOne(tmp, accum, -input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][1], level, param);
        accumulateTrlwe(accum, tmp);
    }
#else
    for (auto i = 0; i < param.n; i++) {
        if (input.a[i] == 0) {
            continue;
        }
        rotateTrlweMinusOne(tmp, accum, input.a[i]);
        externalProductTrgswMPNttInPlace(tmp, bskDft[i][0], level, param);
        accumulateTrlwe(accum, tmp);
    }
#endif
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
        generalExternalProductTrgswMPNtt(accum, trgsws[i], accum, trgsws[i].l, param);
        addTrlev(tmp, accum, temp); // res += input
        accum = std::move(tmp);
    }
}