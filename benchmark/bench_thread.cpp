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
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey{param};
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);


    BootstrappingKeyMPPreRot bskPre{param};
    genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);

    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorus32(pt, param.torusBase);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweFromTorus32(sTlwe, input);
    Trlwe acc{param.k, param.N};
    genNoiselessTrlweSample(acc, v, sTlwe);

    Trlwe out{param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    Tlwe output {param.n};
    vector trgswMPDft(param.n-1, TrgswMPDft{param});

    // rot
    for (auto i = 10; i <= 200; i = i + 10) {
        BENCH500("blindRotateWithPreRotNttMT, bench size=" + to_string(i),
                 blindRotateWithPreRotNttMT(out, trgswMPDft, bskPre.bskFirst, bskPre.bskDft, sTlwe, i, param);)
    }

    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
