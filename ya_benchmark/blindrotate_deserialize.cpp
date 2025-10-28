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
        BootstrappingKeyMP bskMPServer;
        COUNT_TIME("GINX read key", deserializeBskMP(bskMPServer, "bsk_serialized_GINX.bin", param.n);)
        COUNT_TIME("blindRotateGINXNtt", blindRotateJP22Ntt(acc, bskMPServer.bskDft, sTlwe, param);)
    }

    // optimized GINX server procedure
    {
        BootstrappingKeyMPOpt bskMPOptServer;
        COUNT_TIME("GINX_opt read key", deserializeBskMPOpt(bskMPOptServer, "bsk_serialized_GINX_opt_bsk.bin", param.n))
        COUNT_TIME("blindRotateOptNtt", blindRotateOptNtt(out, bskMPOptServer.bskFirst, bskMPOptServer.bskDft, sTlwe, v, param);)
    }

//    // optimized lazy key initialization server procedure
//    {
//        BootstrappingKeyMPLazy bskMPLazyServer;
//        COUNT_TIME("LAZY_PIPE read key", deserializeBskLazyPipe(bskMPLazyServer, "bsk_serialized_PIPE.bin", param.n);)
//        COUNT_TIME("blindRotateLazyPipeNtt", blindRotateLazyPipeNtt(out3, bskMPLazyServer.bskFirst, bskMPLazyServer.bskTrim, bskMPLazyServer.bskDecompA, sTlwe, v, s2Dft, one, param);)
//    }
//
//    // optimized lazy key initialization MT server procedure
//    {
//        BootstrappingKeyMPLazy bskLazy;
//        COUNT_TIME("LAZY_MT read key", deserializeBskLazy(bskLazy, "bsk_serialized_LAZY_MT.bin", param.n);)
//        COUNT_TIME("blindRotateLazyMTNtt", blindRotateLazyMTNtt(out4, bskLazy.bskFirst, bskLazy.bskTrim, bskLazy.bskDecompA, sTlwe, v, s2Dft, param);)
//    }

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
