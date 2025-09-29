#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/initializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    param.l = 2;
    param.batchSize = 40;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    // BootstrappingKey bsKey{param};
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    BootstrappingKeyMP bskMP{param};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMPPreRot bskPre{param};
    genBootstrappingKeyMPPreRotTernary(bskPre, trgswKey, tlweKey, v, param.batchSize, param);


    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorus32(pt, param.torusBase);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey); // encryption
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweFromTorus32(sTlwe, input); // rescale input
    Trlwe acc{param.k, param.N};
    genNoiselessTrlweSample(acc, v, sTlwe); // test polynomial gen
    Trlwe out{param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    Tlwe output {param.n};

    // rot
    BENCH500("blindRotateGINXNtt single thread", blindRotateMPNtt(acc, bskMP.bskDft, sTlwe, param);)
    BENCH500("blindRotateWithPreRotNtt single thread", blindRotateWithPreRotNtt(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);)
    BENCH500("blindRotateGINXNtt multiple threads", blindRotateMPNttMT(acc, bskMP.bskDft, sTlwe, param);)
    BENCH500("blindRotateWithPreRotNtt multiple threads", blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);)

    extractTlweFromTrlwe(tmp, out, param.driftPhase); // sample extract
    switchKeyForTlwe(output, ksKey, tmp, param); // key switch
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase); // decryption
    cout << "decAft: "<< decAft << endl;
    cout << "Correctness: " << ((decAft == pt) ? " Correct (dec = in)." : "Incorrect (dec =/= in).") << endl;
    return 0;
}
