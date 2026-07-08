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
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

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
    generateTestPolynomial(v, param.torusBase, 2 * param.N);


    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    COUNT_TIME("genBootstrappingKeyMPLazy", genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);)

    BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
    symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trlweKey, trlweKey.s, true, param);


    // data gen
    Integer pt = 1;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe acc{param};
    genNoiselessTrlweSample(acc, v, sTlwe);

    Trlwe out2{param};
    Trlwe out31{param};

    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};
    TrgswMPDft one{param};
    encryptTrgswMPNtt(one, 1, trgswKey, 0, param);


    // server side
    // pipelined lazy key initialization with serialization server procedure
    {
        blindRotateLazyPipeSerializationNtt(out2, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft,
            bskMPLazyPipe.bskDecompA, sTlwe, v, bskMPLazyPipeAlt.s2Dft, one,
            "BSK_PIPE.bin", true, param);
    }

    // pipelined lazy key initialization server procedure
    {
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        COUNT_TIME("PIPE_LAZY_INIT blindRotate",
                   blindRotatePipeInitNtt(out31, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft,
                                              bskMPLazyPipeAlt.s2Dft, one, sTlwe, v, "BSK_PIPE.bin", param);)
    }


    // client side
    extractTlweFromTrlwe(tmp, out2, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_SERIA): "<< decAft << endl;
    cout << "err(LAZY_SERIA):" << calTlweError(output, tlweKey, pt) << endl;

    extractTlweFromTrlwe(tmp, out31, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft(LAZY_Pipe_INIT): "<< decAft << endl;
    cout << "err(LAZY_Pipe_INIT):" << calTlweError(output, tlweKey, pt) << endl;
    return 0;
}
