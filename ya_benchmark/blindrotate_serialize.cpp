#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
//    param.N = 1024;
    param.n=4;
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

//    BootstrappingKeyMP bskMP{param, param.lApprox};
//    COUNT_TIME("genBootstrappingKeyMP", genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);)
//    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
//    COUNT_TIME("genBootstrappingKeyMPOpt", genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);)
//    BootstrappingKeyMPOpt bskMPLazy{param, param.lApprox, true};
//    genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);
//    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
//    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);)
//    BootstrappingKeyMPLazy bskMPLazyOpt{param, param.lApprox, true, true};
//    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
    COUNT_TIME("genBootstrappingKeyMPLazyPipeAlt", genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, v, param);)


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
    Trlwe out5{param.k, param.N};
    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};
    TrgswMPDft one{param};
    encryptTrgswMPNtt(one, 1, trgswKey, 0, param);


/*    // server side
    // GINX server procedure
    {
        clearFileCache();
        BootstrappingKeyMP bskMPServer;
        COUNT_TIME("GINX write key", serializeBskMP(bskMP, "BSK_GINX.bin");)
        clearFileCache();
        COUNT_TIME("GINX read key", deserializeBskMP(bskMPServer, "BSK_GINX.bin", param.n);)
        COUNT_TIME("GINX blindRotate", blindRotateJP22Ntt(acc, bskMPServer.bskDft, sTlwe, param);)
    }

    // optimized GINX server procedure
    {
        clearFileCache();
        BootstrappingKeyMPOpt bskMPOptServer;
        COUNT_TIME("GINX_OPT write key", serializeBskMPOpt(bskMPOpt, "BSK_GINX_OPT.bin");)
        clearFileCache();
        COUNT_TIME("GINX_OPT read key", deserializeBskMPOpt(bskMPOptServer, "BSK_GINX_OPT.bin", param.n))
        COUNT_TIME("GINX_OPT blindRotate", blindRotateOptNtt(out, bskMPOptServer.bskFirst, bskMPOptServer.bskDft,
                                                          sTlwe, v, param);)
    }

//    // pipelined lazy key initialization with serialization server procedure
    {
        if (!bskMPLazyPipe.initialized) {
            clearFileCache();
            COUNT_TIME("PIPE_LAZY blindRotate + write key",
                       blindRotateLazyPipeSerializationNtt(out2, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft,
                                                           bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one,
                                                           "BSK_PIPE.bin", true, param);)
            bskMPLazyPipe.initialized = true;
            bskMPLazyPipe.bskDecompA.clear();
        } else {
            blindRotateOptNtt(out2, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft, sTlwe, v, param);
        }
    }*/

//    // pipelined lazy key initialization server procedure
    // {
    //     BootstrappingKeyMPLazyPipe bskMPLazyServer;
    //     if (!bskMPLazyServer.initialized) {
    //         clearFileCache();
    //         COUNT_TIME("PIPE_LAZY read key", deserializeBskLazyPipe(bskMPLazyServer, "BSK_PIPE.bin", param.n);)
    //         COUNT_TIME("PIPE_LAZY blindRotate",
    //                    blindRotateLazyPipeNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft,
    //                                           bskMPLazyServer.bskDecompA, sTlwe, v, s2Dft, one, param);)
    //         bskMPLazyServer.initialized = true;
    //         bskMPLazyServer.bskDecompA.clear();
    //     } else {
    //         blindRotateOptNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, sTlwe, v, param);
    //     }
    // }


/*    // parallel lazy key initialization server procedure
    {
        BootstrappingKeyMPLazy bskMPLazyServer;
        if (!bskMPLazyServer.initialized) {
            clearFileCache();
            COUNT_TIME("PARALLEL_LAZY write key", serializeBskMPLazy(bskMPLazyOpt, "BSK_PAL_LAZY.bin");)
            clearFileCache();
            COUNT_TIME("PARALLEL_LAZY read key", deserializeBskMPLazy(bskMPLazyServer, "BSK_PAL_LAZY.bin", param.n);)
            COUNT_TIME("PARALLEL_LAZY blindRotate",
                       blindRotateLazyMTNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskTrim,
                                            bskMPLazyServer.bskDecompA, sTlwe, v, s2Dft, param);)

            bskMPLazyServer.initialized = true;
            bskMPLazyServer.bskDecompA.clear();
        } else {
            blindRotateOptNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskTrim, sTlwe, v, param);
        }
    }*/

    // pipelined lazy key initialization alternative server procedure
    {
//        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        // clearFileCache();
//            COUNT_TIME("PIPE_LAZY read key", deserializeBskLazyPipe(bskMPLazyServer, "BSK_PIPE.bin", param.n);)
        COUNT_TIME("PIPE_LAZY_ALT blindRotate",
                   blindRotateLazyPipeAltNtt(out5, bskMPLazyPipeAlt.bskFirst, bskMPLazyPipeAlt.bskSecond,
                                             bskMPLazyPipeAlt.bskPrime, sTlwe, v, s2Dft, one, param);)
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
    cout << "decAft(LAZY_SERIA): "<< decAft << endl;
    cout << "err(LAZY_SERIA):" << calTlweError(output, tlweKey, mu) << endl;

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

    extractTlweFromTrlwe(tmp, out5, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_PIPE_ALT): "<< decAft << endl;
    cout << "err(LAZY_PIPE_ALT):" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
