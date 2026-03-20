#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"

int main(int argc, char **argv) {
    YatfheParameters param {};
    param.group = 1;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    COUNT_TIME("genTlweKey", genTlweKey(tlweKey);)
    COUNT_TIME("genTrlweKey", genTrlweKey(trlweKey);)
    COUNT_TIME("genBootstrappingKey", genBootstrappingKey(bsKey, trgswKey, tlweKey, param);)
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    COUNT_TIME("genTlweKeySwitchingKey", genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);)

    Integer plain = 3;
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorusGeneral(mu, param.torusBase, LWE_Q) << endl;
    auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

//    COUNT_TIME("trgswFunctionalBootstrapping", trgswFunctionalBootstrappingNtt(output, input, bsKey, ksKey, v, param);)
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Trlwe accumScaled {param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    COUNT_TIME("rescaleTlweToNewMod", rescaleTlweToNewMod(inputModN2, input);) // rescale to mod 2N
    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(accum, v, inputModN2);) // accum = (X^-b) * (0,...,0,v)
    COUNT_TIME("blindRotateNtt", blindRotateNtt(accum, bsKey.bskDft, inputModN2, param);)
//    COUNT_TIME("blindRotate", blindRotate(accum, bsKey.bsk, inputModN2, param);)
    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, accum, param.driftPhase);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    COUNT_TIME("switchKeyForTlwe", switchKeyForTlwe(output, ksKey, tmp, param);)

    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;

    return 0;
}
