//
// Created by Xintong Song on 2023/12/25.
//
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/key_patterns.h"
#include "yautil/multi_threading.h"

void functionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk.bsk, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotateNtt(accum, bsk.bskDft, inputModN2, param);
    extractTlweFromTrlwe(tmp, accum, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingCrt(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksk.nCurrKey};
    rescaleTlweFromTorus32(inputModN2, input);// rescale to mod 2N
    genNoiselessTrlweSample(tv, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    decompTrlweMcrt(accCRT, tv, param);
    blindRotateApproxCRTNtt(accCRT, bskCRT.bskCRT, inputModN2, param);
    trlweMcrtToCrt(accCRT, param);
    recompTrlweCrt(acc, accCRT, param);
    extractTlweFromTrlwe(tmp, acc, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void genBootstrappingKeyGroup(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    int j = 0;
    const auto group = param.group;
    const auto batchSize = 1 << group;
    const int max = param.n - param.n % param.group;
    for (int i = 0; i < max; i = i + group) {
        int combined = 0;
        for (int i2 = 0; i2 < group; i2++) {
            const auto s = tlweKey.s[i + i2] << (group - 1 - i2);
            combined |= s;
        }

        for (int k = 0; k < batchSize; ++k) {
            const int mu = combined == k ? 1 : 0;
            encryptTrgswNtt(bsk.bsk[j + k], bsk.bskDft[j + k], mu, trgswKey, 0, param);
        }
        j += batchSize;
    }
    for (int i = max; i < tlweKey.n; ++i) {
        encryptTrgswNtt(bsk.bsk[j], bsk.bskDft[j], tlweKey.s[i], trgswKey, 0, param);
        j++;
    }
}

void genBootstrappingKeyNormal(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    for (auto i = 0; i < bsk.n; i++) {
        encryptTrgswNtt(bsk.bsk[i], bsk.bskDft[i], tlweKey.s[i], trgswKey, 0, param);
    }
}

void genBootstrappingKey(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    if (bsk.group == 1) {
        genBootstrappingKeyNormal(bsk, trgswKey, tlweKey, param);
        return;
    }
    genBootstrappingKeyGroup(bsk, trgswKey, tlweKey, param);
}

void genBootstrappingKeyMP(BootstrappingKeyMP& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    for (auto i = 0; i < bsk.n; i++) {
#ifdef TERNARY
        const auto si = tlweKey.s[i];
        if(si == 0) {
            encryptTrgswMPNtt(bsk.bskDft[i][0], 0, trgswKey, 0, param);
            encryptTrgswMPNtt(bsk.bskDft[i][1], 0, trgswKey, 0, param);
        } else if (si == 1) {
            encryptTrgswMPNtt(bsk.bskDft[i][0], 1, trgswKey, 0, param);
            encryptTrgswMPNtt(bsk.bskDft[i][1], 0, trgswKey, 0, param);
        } else {
            encryptTrgswMPNtt(bsk.bskDft[i][0], 0, trgswKey, 0, param);
            encryptTrgswMPNtt(bsk.bskDft[i][1], 1, trgswKey, 0, param);
        }
#else
        encryptTrgswMPNtt(bsk.bskDft[i][0], tlweKey.s[i], trgswKey, 0, param);
#endif
    }
}

void genBootstrappingKeyMPPreRot(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                 const TorusPolynomial& v, const int batchSize, const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n; start += batchSize) {
        futures.clear();
        const auto minVal = min(start + batchSize, n);
        for (int i = start; i < minVal; i++) {
            if (i == 0) {
                if (tlweKey.s[i] == 1) {
                    symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
                    symEncTrlweSingleSample(bsk.bskFirst[1], trgswKey.trlweKey, 0);
                } else {
                    symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
                    symEncTrlweMultiSample(bsk.bskFirst[1], trgswKey.trlweKey, v.coeffs);
                }
                continue;
            }
            futures.emplace_back(pool.enqueue([i, &tlweKey, &bsk, &trgswKey, &param] {
                auto j = i-1;
                if (tlweKey.s[i] == 1) {
                    encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                } else {
                    encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
                }
            }
            ));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

void genBootstrappingKeyMPPreRotTernary(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                        const TorusPolynomial& v, const int batchSize, const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n; start += batchSize) {
        futures.clear();
        const auto minVal = min(start + batchSize, n);
        for (int i = start; i < minVal; i++) {
            if (i == 0) {
                if (tlweKey.s[i] == 1) {
                    symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
                    symEncTrlweSingleSample(bsk.bskFirst[1], trgswKey.trlweKey, 0);
                    symEncTrlweSingleSample(bsk.bskFirst[2], trgswKey.trlweKey, 0);
                } else if (tlweKey.s[i] == 0) {
                    symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
                    symEncTrlweSingleSample(bsk.bskFirst[1], trgswKey.trlweKey, 0);
                    symEncTrlweMultiSample(bsk.bskFirst[2], trgswKey.trlweKey, v.coeffs);
                } else {
                    symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
                    symEncTrlweMultiSample(bsk.bskFirst[1], trgswKey.trlweKey, v.coeffs);
                    symEncTrlweSingleSample(bsk.bskFirst[2], trgswKey.trlweKey, 0);
                }
                continue;
            }
            futures.emplace_back(pool.enqueue([i, &tlweKey, &bsk, &trgswKey, &param] {
                auto j = i-1;
                if (tlweKey.s[i] == 1) {
                    encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][2],0, trgswKey, 0, param);
                } else if (tlweKey.s[i] == 0) {
                    encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][2],1, trgswKey, 0, param);
                } else {
                    encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
                    encryptTrgswMPNtt(bsk.bskDft[j][2],0, trgswKey, 0, param);
                }
            }
            ));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

void genBootstrappingKeyInternal(BootstrappingKeyInternal& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                 const YatfheParameters& param) {
    const auto n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n); i++) {
            const auto j = i / 2;
            futures.emplace_back(pool.enqueue([i, j, &tlweKey, &bsk, &trgswKey, &param] {
                if (i % 2 == 0) {
                    if (tlweKey.s[i] == 1) {
                        encryptTrgswMP(bsk.bsk[j][0],1, trgswKey, 0, param);
                        encryptTrgswMP(bsk.bsk[j][1],0, trgswKey, 0, param);
                    } else {
                        encryptTrgswMP(bsk.bsk[j][0],0, trgswKey, 0, param);
                        encryptTrgswMP(bsk.bsk[j][1],1, trgswKey, 0, param);
                    }
                } else {
                    if (tlweKey.s[i] == 1) {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                    } else {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
                    }
                }
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
    // for (auto i = 0; i < param.n; i++) {
    //     const auto j = i / 2;
    //     if (i % 2 == 0) {
    //         if (tlweKey.s[i] == 1) {
    //             encryptTrgswMP(bsk.bsk[j][0],1, trgswKey, 0, param);
    //             encryptTrgswMP(bsk.bsk[j][1],0, trgswKey, 0, param);
    //         } else {
    //             encryptTrgswMP(bsk.bsk[j][0],0, trgswKey, 0, param);
    //             encryptTrgswMP(bsk.bsk[j][1],1, trgswKey, 0, param);
    //         }
    //     } else {
    //         if (tlweKey.s[i] == 1) {
    //             encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
    //             encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
    //         } else {
    //             encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
    //             encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
    //         }
    //     }
    // }
}

void genBootstrappingKeyInternalAsym(BootstrappingKeyInternalAsym& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                     const YatfheParameters& param) {
    const auto n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n); i++) {
            const auto j = i / 2;
            futures.emplace_back(pool.enqueue([i, j, &tlweKey, &bsk, &trgswKey, &param] {
                if (i % 2 != 0) {
                    TorusPolynomial s{param.N};
                    s.coeffs[0] = 1;
                    if (tlweKey.s[i] == 1) {
                        encTrlevMultiSample(bsk.bsk[j][0], trgswKey.trlweKey, s, param);
                        encTrlevSingleSample(bsk.bsk[j][1], trgswKey.trlweKey, 0, param);
                    } else {
                        encTrlevSingleSample(bsk.bsk[j][0], trgswKey.trlweKey, 0, param);
                        encTrlevMultiSample(bsk.bsk[j][1], trgswKey.trlweKey, s, param);
                    }
                } else {
                    if (tlweKey.s[i] == 1) {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                    } else {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
                    }
                }
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

void genBootstrappingKeyInternalAsymOpt(BootstrappingKeyInternalAsymOpt& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                        const TorusPolynomial& v, const YatfheParameters& param) {
    const auto n = param.n;
    const int batchSize = 32;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
    for (int start = 0; start < n; start += batchSize) {
        futures.clear();
        for (int i = start; i < min(start + batchSize, n); i++) {
            if (i == n - 1) {
                if (tlweKey.s[i] == 1) {
                    symEncTrlweMultiSample(bsk.bskLast[0], trgswKey.trlweKey, v.coeffs);
                    symEncTrlweSingleSample(bsk.bskLast[1], trgswKey.trlweKey, 0);
                } else {
                    symEncTrlweSingleSample(bsk.bskLast[0], trgswKey.trlweKey, 0);
                    symEncTrlweMultiSample(bsk.bskLast[1], trgswKey.trlweKey, v.coeffs);
                }
                continue;
            }
            const auto j = i / 2;
            futures.emplace_back(pool.enqueue([i, j, &tlweKey, &bsk, &trgswKey, &param] {
                if (i % 2 != 0) {
                    TorusPolynomial s{param.N};
                    s.coeffs[0] = 1;
                    if (tlweKey.s[i] == 1) {
                        encTrlevMultiSample(bsk.bsk[j][0], trgswKey.trlweKey, s, param);
                        encTrlevSingleSample(bsk.bsk[j][1], trgswKey.trlweKey, 0, param);
                    } else {
                        encTrlevSingleSample(bsk.bsk[j][0], trgswKey.trlweKey, 0, param);
                        encTrlevMultiSample(bsk.bsk[j][1], trgswKey.trlweKey, s, param);
                    }
                } else {
                    if (tlweKey.s[i] == 1) {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],1, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],0, trgswKey, 0, param);
                    } else {
                        encryptTrgswMPNtt(bsk.bskDft[j][0],0, trgswKey, 0, param);
                        encryptTrgswMPNtt(bsk.bskDft[j][1],1, trgswKey, 0, param);
                    }
                }
            }));
        }
        for (auto& f : futures) {
            f.wait();
        }
    }
}

void genBootstrappingKeyApproxCrt(BootstrappingKeyCRT& bskCRT, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param) {
    BootstrappingKey bsk{param};
    for (auto i = 0; i < bsk.n; i++) {
        encryptTrgswApproxCRT(bsk.bsk[i], param, trgswKey, tlweKey.s[i]);
    }
    decompBootstrappingKeyMcrt(bskCRT, bsk, param);
}

void decompBootstrappingKeyMcrt(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param) {
    for (size_t i = 0; i < param.n; i++) {
        decompTrgswMcrt(bskCRT.bsk8[i], bsk.bsk[i], param);
        auto& bsk8D = bskCRT.bsk8[i];
        auto& bskNttD = bskCRT.bskCRT[i];
        for (size_t d = 0; d < param.d; d++) {
            auto& bsk8 = bsk8D[d];
            auto& bskNtt = bskNttD[d];
            NttNative24::applyNttForRgsw(bskNtt, bsk8);
        }
    }
}