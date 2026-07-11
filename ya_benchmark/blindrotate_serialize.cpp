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
    param.batchSize = 3;
    param.tasksPerThread = 1;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, T:%d, b:%d, l:%d, lA:%d\n", param.n, param.k, param.N, param.torusBits, param.radixBits, param.l, param.lApprox);

    // client side
    // key gen
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.errorB = param.rlweNoiseB;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TorusPolynomial v {param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMP bskMP{param, param.lApprox};
    COUNT_TIME("genBootstrappingKeyMP", genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);)
    BootstrappingKeyWWL24 bskWWL24{param, param.lApprox};
    // BootstrappingKeyWWL24 bskWWL24{param, param.lApprox + 1};
    COUNT_TIME("genBootstrappingKeyWWL24", genBootstrappingKeyWWL24(bskWWL24, trgswKey, tlweKey, param);)
    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    COUNT_TIME("genBootstrappingKeyMPOpt", genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPOpt bskMPLazy{param, param.lApprox, true};
    genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPLazy bskMPLazyOpt{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);)
    BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
    symEncTrlevWithKeyNtt(bskMPLazyPipe.s2Dft, trlweKey, trlweKey.s, true, param);
    symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trlweKey, trlweKey.s, true, param);
    COUNT_TIME("genBootstrappingKeyMPLazyPipeAlt", genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, v, param);)


    // data gen
    Integer pt = 0;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);

    Trlwe acc{param};
    Trlwe out{param};
    Trlwe out2{param};
    Trlwe out3{param};
    Trlwe out31{param};
    Trlwe out4{param};
    Trlwe out5{param};
    Trlwe out6{param};
    Trlwe out7{param};
    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};


    // server side
    // GINX server procedure
    {
        clearFileCache();
        BootstrappingKeyMP bskMPServer;
        COUNT_TIME("GINX write key", serializeBskMP(bskMP, "BSK_GINX.bin");)
        clearFileCache();
        COUNT_TIME("GINX read key", deserializeBskMP(bskMPServer, "BSK_GINX.bin", param.n);)
        genNoiselessTrlweSample(acc, v, sTlwe);
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


    // pipelined lazy key initialization server procedure
    {
        clearFileCache();
        COUNT_TIME("PIPE_LAZY_INIT write key", serializeBskLazyPipe(bskMPLazyPipe, "BSK_PIPE_INIT.bin");)
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        clearFileCache();
        COUNT_TIME("PIPE_LAZY_INIT blindRotate",
                   blindRotatePipeInitNtt(out31, bskMPLazyServer, sTlwe, v, "BSK_PIPE_INIT.bin", param);)
    }

    // pipelined lazy key initialization server procedure
    {
        clearFileCache();
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        COUNT_TIME("deserializeBskLazyPipe", deserializeBskLazyPipe(bskMPLazyServer, "BSK_PIPE_INIT.bin", param.n);)
        clearFileCache();
        COUNT_TIME("PIPE_LAZY blindRotate", blindRotateLazyPipeNtt(out3, bskMPLazyServer, sTlwe, v, param);)
    }


    // parallel lazy key initialization server procedure
    {
        BootstrappingKeyMPLazy bskMPLazyServer;
        if (!bskMPLazyServer.initialized) {
            clearFileCache();
            COUNT_TIME("PARALLEL_LAZY write key", serializeBskMPLazy(bskMPLazyOpt, "BSK_PAL_LAZY.bin");)
            clearFileCache();
            COUNT_TIME("PARALLEL_LAZY read key", deserializeBskMPLazy(bskMPLazyServer, "BSK_PAL_LAZY.bin", param.n);)
            COUNT_TIME("PARALLEL_LAZY blindRotate",
                       blindRotateLazyMTNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskTrim,
                                            bskMPLazyServer.bskDecompA, sTlwe, v, bskMPLazyPipeAlt.s2Dft, param);)

            // bskMPLazyServer.initialized = true;
            // bskMPLazyServer.bskDecompA.clear();
        } else {
            blindRotateOptNtt(out4, bskMPLazyServer.bskFirst, bskMPLazyServer.bskTrim, sTlwe, v, param);
        }
    }

    // pipelined lazy key initialization alternative server procedure
    {
        clearFileCache();
        COUNT_TIME("PIPE_LAZY_ALT blindRotate",
                   blindRotateLazyPipeAltNtt(out5, bskMPLazyPipeAlt.bskFirst, bskMPLazyPipeAlt.bskPrime,
                       bskMPLazyPipeAlt.s2Dft, sTlwe, v, param);)
        COUNT_TIME("PIPE_LAZY_ALT write key",serializeBskLazyPipeAlt(bskMPLazyPipeAlt, "BSK_PIPE_ALT.bin");)
    }

    {
        BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAltServer;
        clearFileCache();
        COUNT_TIME("PIPE_LAZY_ALT_INIT blindRotate",
                   blindRotateLazyPipeAltInitNtt(out6, bskMPLazyPipeAltServer.bskFirst, bskMPLazyPipeAltServer.bskPrime,
                       bskMPLazyPipeAltServer.s2Dft, sTlwe, v, "BSK_PIPE_ALT.bin", param);)
    }

    // WWL24 procedure
    {
        clearFileCache();
        BootstrappingKeyWWL24 bskWWL24Server;
        COUNT_TIME("WWL24 write key", serializeBskWWL24(bskWWL24, "BSK_WWL.bin");)
        clearFileCache();
        COUNT_TIME("WWL24 read key", deserializeBskWWL24(bskWWL24Server, "BSK_WWL.bin", param.n);)
        genNoiselessTrlweSample(out7, v, sTlwe);
        COUNT_TIME("WWL24 blindRotate", blindRotateWWL24Ntt(out7, bskWWL24Server.bskDft, sTlwe, bskWWL24Server.s2Dft, param);)
    }

    // client side
    extractTlweFromTrlwe(tmp, acc, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(GINX): "<< decAft << endl;
    cout << "err(GINX):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(GINX_OPT): "<< decAft << endl;
    cout << "err(GINX_OPT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out3, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe): "<< decAft << endl;
    cout << "err(LAZY_Pipe):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out31, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe_INIT): "<< decAft << endl;
    cout << "err(LAZY_Pipe_INIT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out4, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_MT): "<< decAft << endl;
    cout << "err(LAZY_MT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out5, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_PIPE_ALT): "<< decAft << endl;
    cout << "err(LAZY_PIPE_ALT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out6, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(PIPE_LAZY_ALT_INIT): "<< decAft << endl;
    cout << "err(PIPE_LAZY_ALT_INIT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out7, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(WWL24): "<< decAft << endl;
    cout << "err(WWL24):" << calTlweError(output, tlweKey, pt) << endl;

    return 0;
}
