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

    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);)
    symEncTrlevWithKeyNtt(bskMPLazyPipe.s2Dft, trlweKey, trlweKey.s, true, param);


    // data gen
    Integer pt = 0;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);

    Trlwe out{param};
    Trlwe out1{param};
    Trlwe out2{param};
    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};

    {
        genNoiselessTrlweSample(out, v, sTlwe);
        COUNT_TIME("GINX blindRotate", blindRotateJP22Ntt(out, bskMP.bskDft, sTlwe, param);)
    }
    {
        serializeBskLazyPipe(bskMPLazyPipe, "BSK_PIPE.bin");
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        COUNT_TIME("PIPE_LAZY_INIT blindRotate",
                   blindRotatePipeInitNtt(out1, bskMPLazyServer, sTlwe, v, "BSK_PIPE.bin", param);)
    }
    {
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        deserializeBskLazyPipe(bskMPLazyServer, "BSK_PIPE.bin", param.n);
        COUNT_TIME("PIPE_LAZY blindRotate", blindRotateLazyPipeNtt(out2, bskMPLazyServer, sTlwe, v, param);)
    }
    //
    IntPolynomial p{param.N};
    symDecTrlweToInt(p, out, trlweKey, param.torusBase);
    printArray(p.coeffs, "p");


    // client side
    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(GINX): "<< decAft << endl;
    cout << "err(GINX):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out1, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe_INIT): "<< decAft << endl;
    cout << "err(LAZY_Pipe_INIT):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out2, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe): "<< decAft << endl;
    cout << "err(LAZY_Pipe):" << calTlweError(output, tlweKey, pt) << endl;

    return 0;
}
