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
    initYatfhe(param);

    // key gen
    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    BootstrappingKeyInternal bsk {param};
    genBootstrappingKeyInternal(bsk, trgswKey, tlweKey, param);
    BootstrappingKey bskNor {param};
    genBootstrappingKey(bskNor, trgswKey, tlweKey, param);
    BootstrappingKeyInternalAsym bskAsym {param};
    genBootstrappingKeyInternalAsym(bskAsym, trgswKey, tlweKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyInternalAsymOpt bskAsymOpt {param};
    genBootstrappingKeyInternalAsymOpt(bskAsymOpt, trgswKey, tlweKey, v, param);

    TrlevDft s2Dft(param);
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

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

    // rot
    COUNT_TIME("blindRotateNtt", blindRotateNtt(acc, bskNor.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateInternalPireWiseNtt", blindRotateInternalPairWiseNtt(acc, bsk.bsk, bsk.bskDft, sTlwe, param);)
    COUNT_TIME("blindRotateInternalPairWiseAsymNtt", blindRotateInternalPairWiseAsymNtt(acc, bskAsym.bsk, bskAsym.bskDft, sTlwe, s2Dft, param);)
    COUNT_TIME("blindRotateInternalPairWiseAsymOptNtt", blindRotateInternalPairWiseAsymOptNtt(out, bskAsymOpt.bsk, bskAsymOpt.bskLast, bskAsymOpt.bskDft, sTlwe, s2Dft, param);)

    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
