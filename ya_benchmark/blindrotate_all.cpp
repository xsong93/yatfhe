#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    param.N = 1024;
    param.batchSize = 3;
    param.tasksPerThread = 1;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    // client side
    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TrlevDft s2Dft{param, param.l};
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMP bskMP{param, param.lApprox};
    COUNT_TIME("genBootstrappingKeyMP", genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);)
    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    COUNT_TIME("genBootstrappingKeyMPOpt", genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPOpt bskMPLazy{param, param.lApprox, true};
    genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPLazy bskMPLazyOpt{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);)


    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe acc{param.k, param.N};
    genNoiselessTrlweSample(acc, v, sTlwe);

    Trlwe out{param.k, param.N};
    Trlwe out2{param.k, param.N};
    Trlwe out3{param.k, param.N};
    Trlwe out4{param.k, param.N};
    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};
    TrgswMPDft one{param};
    encryptTrgswMPNtt(one, 1, trgswKey, 0, param);


    // server side
    // GINX server procedure
    {
        BootstrappingKeyMP bskMPServer{param, param.lApprox};
        COUNT_TIME("blindRotateGINXNtt key", bskMPServer = bskMP;)
        COUNT_TIME("GINX serialize", {
            vector<char> buffer;
            for (auto i = 0; i < param.n; i++) {
                std::ostringstream oss(std::ios::binary);
                serialize(bskMPServer.bskDft[i][0], oss);
                const std::string &data = oss.str();
                buffer.insert(buffer.end(), data.begin(), data.end());
            }
            std::ofstream outFile("bsk_serialized_GINX.bin", std::ios::binary | std::ios::app);
            outFile.write(buffer.data(), buffer.size());
            buffer.clear();
        })
        COUNT_TIME("blindRotateGINXNtt", blindRotateJP22Ntt(acc, bskMPServer.bskDft, sTlwe, param);)
    }

    // optimized GINX server procedure
    {
        BootstrappingKeyMPOpt bskMPOptServer{param, param.lApprox, false};
        COUNT_TIME("blindRotateOptNtt key", bskMPOptServer = bskMPOpt;)
        COUNT_TIME("GINX_opt serialize", {
            vector<char> buffer;
            std::ostringstream oss1(std::ios::binary);
            serialize(bskMPOptServer.bskFirst[0], oss1);
            const std::string &data1 = oss1.str();
            buffer.insert(buffer.end(), data1.begin(), data1.end());
            for (auto i = 0; i < param.n - 1; i++) {
                std::ostringstream oss(std::ios::binary);
                serialize(bskMPOptServer.bskDft[i][0], oss);
                const std::string &data = oss.str();
                buffer.insert(buffer.end(), data.begin(), data.end());
            }
            std::ofstream outFile2("bsk_serialized_GINX_opt_bsk.bin", std::ios::binary | std::ios::app);
            outFile2.write(buffer.data(), buffer.size());
            buffer.clear();
        })
        COUNT_TIME("blindRotateOptNtt", blindRotateOptNtt(out, bskMPOptServer.bskFirst, bskMPOptServer.bskDft, sTlwe, v, param);)
    }

    // parallel naive lazy key initialization server procedure
    {
        BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true};
        COUNT_TIME("blindRotateLazyNtt key", bskMPLazyServer = bskMPLazy;)
        COUNT_TIME("blindRotateLazyNtt", blindRotateLazyNtt(out2, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyServer.initialized, sTlwe, v, s2Dft, param);)
    }

    // optimized lazy key initialization server procedure
    {
        BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true, true};
        if (!bskMPLazyServer.initialized) {
            COUNT_TIME("blindRotateLazyPipeNtt key", {
                bskMPLazyServer.bskFirst = bskMPLazyPipe.bskFirst;
                bskMPLazyServer.bskDft[0] = bskMPLazyPipe.bskFull[0];
                bskMPLazyServer.bskDft[1] = bskMPLazyPipe.bskFull[1];
                for (auto i = 2; i < param.n - 1; i++) {
                    bskMPLazyServer.bskDft[i] = bskMPLazyPipe.bskTrim[i - 2];
                }
            })
            COUNT_TIME("blindRotateLazyPipeNtt", blindRotateLazyPipeNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one, param);)
            bskMPLazyServer.initialized = true;
        } else {
            blindRotateOptNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, sTlwe, v, param);
        }
    }

    // optimized lazy key initialization with serialization server procedure
    {
        BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true, true};
        if (!bskMPLazyServer.initialized) {
            COUNT_TIME("blindRotateLazyPipeNtt key", {
                bskMPLazyServer.bskFirst = bskMPLazyPipe.bskFirst;
                bskMPLazyServer.bskDft[0] = bskMPLazyPipe.bskFull[0];
                bskMPLazyServer.bskDft[1] = bskMPLazyPipe.bskFull[1];
                for (auto i = 2; i < param.n - 1; i++) {
                    bskMPLazyServer.bskDft[i] = bskMPLazyPipe.bskTrim[i - 2];
                }
            })
            COUNT_TIME("blindRotateLazyPipeSerializationNtt", blindRotateLazyPipeSerializationNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one, param);)
            bskMPLazyServer.initialized = true;
        } else {
            blindRotateOptNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, sTlwe, v, param);
        }
    }


    // optimized lazy key initialization MT server procedure
    {
        BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true, true};
        if (!bskMPLazyServer.initialized) {
            COUNT_TIME("blindRotateLazyMTNtt key", {
                bskMPLazyServer.bskFirst = bskMPLazyOpt.bskFirst;
                for (auto i = 0; i < param.n - 1; i++) {
                    bskMPLazyServer.bskDft[i] = bskMPLazyOpt.bskTrim[i];
                }
            })
            COUNT_TIME("blindRotateLazyMTNtt", blindRotateLazyMTNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyOpt.bskDecompA, sTlwe, v, s2Dft, param);)
            COUNT_TIME("LAZY_MT serialize", {
                vector<char> buffer;
                std::ostringstream oss1(std::ios::binary);
                serialize(bskMPLazyOpt.bskFirst[0], oss1);
                const std::string &data1 = oss1.str();
                buffer.insert(buffer.end(), data1.begin(), data1.end());
                for (auto i = 0; i < param.n - 1; i++) {
                    std::ostringstream oss(std::ios::binary);
                    serialize(bskMPLazyOpt.bskTrim[i][0], oss);
                    const std::string &data = oss.str();
                    buffer.insert(buffer.end(), data.begin(), data.end());
                }
                for (auto i = 0; i < param.n - 1; i++) {
                    std::ostringstream oss(std::ios::binary);
                    serializeNestedVector(bskMPLazyOpt.bskDecompA[i][0], oss);
                    const std::string &data = oss.str();
                    buffer.insert(buffer.end(), data.begin(), data.end());
                }
                std::ofstream outFile2("bsk_serialized_LAZY_MT.bin", std::ios::binary | std::ios::app);
                outFile2.write(buffer.data(), buffer.size());
                buffer.clear();
            })
            bskMPLazyServer.initialized = true;
        } else {
            blindRotateOptNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, sTlwe, v, param);
        }
    }

    // client side
    extractTlweFromTrlwe(tmp, acc, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(GINX): "<< decAft << endl;
    cout << "err(GINX):" << calTlweError(output, tlweKey, mu) << endl;

    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(GINX_OPT): "<< decAft << endl;
    cout << "err(GINX_OPT):" << calTlweError(output, tlweKey, mu) << endl;

    extractTlweFromTrlwe(tmp, out2, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY): "<< decAft << endl;
    cout << "err(LAZY):" << calTlweError(output, tlweKey, mu) << endl;

    extractTlweFromTrlwe(tmp, out3, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe): "<< decAft << endl;
    cout << "err(LAZY_Pipe):" << calTlweError(output, tlweKey, mu) << endl;

    extractTlweFromTrlwe(tmp, out4, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_MT): "<< decAft << endl;
    cout << "err(LAZY_MT):" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
