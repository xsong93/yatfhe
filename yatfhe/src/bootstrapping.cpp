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
    Trlwe accumScaled {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweToNewMod(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotate(accum, bsk.bsk, inputModN2, param);
    rescaleTrlweToNewMod(accumScaled, accum, LWE_Q, TORUS_Q);
    extractTlweFromTrlwe(tmp, accumScaled, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Trlwe accumScaled {param.k, param.N};
    Tlwe tmp {ksk.nCurrKey};
    rescaleTlweToNewMod(inputModN2, input); // rescale to mod 2N
    genNoiselessTrlweSample(accum, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    blindRotateNtt(accum, bsk.bskDft, inputModN2, param);
    rescaleTrlweToNewMod(accumScaled, accum, LWE_Q, TORUS_Q);
    extractTlweFromTrlwe(tmp, accumScaled, param.driftPhase); // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    switchKeyForTlwe(out, ksk, tmp, param);
}

void functionalBootstrappingCrt(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param) {
    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    Trlwe accumScaled {param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksk.nCurrKey};
    rescaleTlweToNewMod(inputModN2, input);// rescale to mod 2N
    genNoiselessTrlweSample(tv, v, inputModN2); // accum = (X^-b) * (0,...,0,v)
    decompTrlweMcrt(accCRT, tv, param);
    blindRotateApproxCRTNtt(accCRT, bskCRT.bskCRT, inputModN2, param);
    trlweMcrtToCrt(accCRT, param);
    recompTrlweCrt(acc, accCRT, param);
    rescaleTrlweToNewMod(accumScaled, acc, LWE_Q, TORUS_Q);
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

void genBootstrappingKeyMPOpt(BootstrappingKeyMPOpt& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const TorusPolynomial& v, const YatfheParameters& param) {
    {
#ifdef TERNARY
        if (tlweKey.s[0] == 0) {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
            symEncTrlweSingleSample(bsk.bskFirst[1], trgswKey.trlweKey, 0);
        } else if (tlweKey.s[0] == 1) {
            symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
            symEncTrlweSingleSample(bsk.bskFirst[1], trgswKey.trlweKey, 0);
        } else {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
            symEncTrlweMultiSample(bsk.bskFirst[1], trgswKey.trlweKey, v.coeffs);
        }
#else
        if (tlweKey.s[0] == 1) {
            symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
        } else {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
        }
#endif
    }

    for (auto i = 0; i < bsk.n - 1; i++) {
#ifdef TERNARY
        const auto si = tlweKey.s[i+1];
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
        encryptTrgswMPNtt(bsk.bskDft[i][0], tlweKey.s[i+1], trgswKey, 0, param);
#endif
    }
}

void genBootstrappingKeyMPLazy(BootstrappingKeyMPLazy& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                               const TorusPolynomial& v, const YatfheParameters& param) {
    // process first key components
    {
#ifdef TERNARY
#else
        if (tlweKey.s[0] == 1) {
            symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
        } else {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
        }
#endif
    }

    // process remaining n - 1 components
    for (auto i = 0; i < bsk.n - 1; i++) {
#ifdef TERNARY
#else
        TrgswMP tmp{param, bsk.level, true};
        encryptTrgswMP(tmp, tlweKey.s[i + 1], trgswKey, 0, param);
        for (auto l0 = 0; l0 < bsk.level; l0++) {
            auto& a = tmp.cPrime[l0].a;
            auto& dA = bsk.bskDecompA[i][0][l0];
            for (auto k = 0; k < param.k; k++) {
                for (auto j = 0; j < param.N; j++) {
                    DecomposedData d{param.l};
                    gadgetDecompose(d, a[k].coeffs[j], param);
                    for (auto l = 0; l < param.l; l++) {
                        dA[l][k].coeffs[j] = d.value[l] * d.sign;
                    }
                }
            }
            NttHexl::applyNtt(bsk.bskTrim[i][0].cPrime[l0].b, tmp.cPrime[l0].b);
        }
#endif
    }
}

void genBootstrappingKeyMPLazyPipe(BootstrappingKeyMPLazyPipe& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                                   const TorusPolynomial& v, const YatfheParameters& param) {
    // process first two key components
    {
#ifdef TERNARY
#else
        if (tlweKey.s[0] == 1) {
            symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
        } else {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
        }
        encryptTrgswMPNtt(bsk.bskDft[0][0], tlweKey.s[1], trgswKey, 0, param);
        encryptTrgswMPNtt(bsk.bskDft[1][0], tlweKey.s[2], trgswKey, 0, param);
#endif
    }

    // process remaining n - 3 components
    for (auto i = 0; i < bsk.n - 3; i++) {
#ifdef TERNARY
#else
        TrgswMP tmp{param, bsk.level, true};
        encryptTrgswMP(tmp, tlweKey.s[i + 3], trgswKey, 0, param);
        for (auto l0 = 0; l0 < bsk.level; l0++) {
            auto& a = tmp.cPrime[l0].a;
            auto& dA = bsk.bskDecompA[i][0][l0];
            for (auto k = 0; k < param.k; k++) {
                for (auto j = 0; j < param.N; j++) {
                    DecomposedData d{param.l};
                    gadgetDecompose(d, a[k].coeffs[j], param);
                    for (auto l = 0; l < param.l; l++) {
                        dA[l][k].coeffs[j] = d.value[l] * d.sign;
                    }
                }
            }
            NttHexl::applyNtt(bsk.bskDft[i + 2][0].cPrime[l0].b, tmp.cPrime[l0].b);
        }
#endif
    }
}

void genBootstrappingKeyMPLazyPipeAlt(BootstrappingKeyMPLazyPipeAlt& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                   const TorusPolynomial& v, const YatfheParameters& param) {
    // process first two key components
    {
#ifdef TERNARY
#else
        if (tlweKey.s[0] == 1) {
            symEncTrlweMultiSample(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs);
        } else {
            symEncTrlweSingleSample(bsk.bskFirst[0], trgswKey.trlweKey, 0);
        }
        // encryptTrgswMPNtt(bsk.bskSecond[0], tlweKey.s[1], trgswKey, 0, param);
#endif
    }

    // process remaining n - 2 components
    for (auto i = 0; i < bsk.n - 1; i++) {
#ifdef TERNARY
#else
        encryptTrgswMP(bsk.bskPrime[i][0], tlweKey.s[i + 1], trgswKey, 0, param);
#endif
    }
}

void genBootstrappingKeyMPFixedNoise(BootstrappingKeyMP& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const Torus noise, const YatfheParameters& param) {
    for (auto i = 0; i < bsk.n; i++) {
#ifdef TERNARY
        const auto si = tlweKey.s[i];
        if(si == 0) {
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][0], 0, trgswKey, 0, noise, param);
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][1], 0, trgswKey, 0, noise, param);
        } else if (si == 1) {
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][0], 1, trgswKey, 0, noise, param);
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][1], 0, trgswKey, 0, noise, param);
        } else {
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][0], 0, trgswKey, 0,noise, param);
            encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][1], 1, trgswKey, 0, noise, param);
        }
#else
        encryptTrgswMPFixedNoiseNtt(bsk.bskDft[i][0], tlweKey.s[i], trgswKey, 0, noise, param);
#endif
    }
}

void genBootstrappingKeyMPPreRot(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                 const TorusPolynomial& v, const int batchSize, const YatfheParameters& param) {
    const auto n = param.n;
    auto& pool = ThreadPool::instance();
    vector<future<void>> futures;
    futures.reserve(batchSize);
#ifdef TERNARY
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
#else
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
#endif
}

void genBootstrappingKeyMPPreRotTernaryFixedNoise(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                                  const TorusPolynomial& v, const int batchSize, const Torus noise, const YatfheParameters& param) {
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
                    symEncTrlweMultiSampleFixedNoise(bsk.bskFirst[0], trgswKey.trlweKey, v.coeffs, noise);
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[1], trgswKey.trlweKey, 0, noise);
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[2], trgswKey.trlweKey, 0, noise);
                } else if (tlweKey.s[i] == 0) {
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[0], trgswKey.trlweKey, 0, noise);
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[1], trgswKey.trlweKey, 0, noise);
                    symEncTrlweMultiSampleFixedNoise(bsk.bskFirst[2], trgswKey.trlweKey, v.coeffs, noise);
                } else {
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[0], trgswKey.trlweKey, 0, noise);
                    symEncTrlweMultiSampleFixedNoise(bsk.bskFirst[1], trgswKey.trlweKey, v.coeffs, noise);
                    symEncTrlweSingleSampleFixedNoise(bsk.bskFirst[2], trgswKey.trlweKey, 0, noise);
                }
                continue;
            }
            futures.emplace_back(pool.enqueue([i, noise, &tlweKey, &bsk, &trgswKey, &param] {
                auto j = i-1;
                if (tlweKey.s[i] == 1) {
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][0],1, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][1],0, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][2],0, trgswKey, 0, noise, param);
                } else if (tlweKey.s[i] == 0) {
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][0],0, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][1],0, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][2],1, trgswKey, 0, noise, param);
                } else {
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][0],0, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][1],1, trgswKey, 0, noise, param);
                    encryptTrgswMPFixedNoiseNtt(bsk.bskDft[j][2],0, trgswKey, 0, noise, param);
                }
            }
            ));
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