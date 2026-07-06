#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    param.q = Q_CRT;
    param.torusBits = 32;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKeyCRT bsKeyCRT{param};
    TlweKeySwitchingKey ksKey{param};
    COUNT_TIME("genTlweKey", genTlweKey(tlweKey);)
    COUNT_TIME("genTrlweKey", genTrlweKey(trlweKey);)
    COUNT_TIME("genBootstrappingKeyApproxCrt", genBootstrappingKeyApproxCrt(bsKeyCRT, trgswKey, tlweKey, param);)
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    COUNT_TIME("genTlweKeySwitchingKey", genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);)

    Integer plain = 3;
    Torus mu = modSwitchToTorusGeneral(plain, param.torusBase, LWE_Q);
    TorusPolynomial v{param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    Tlwe input{param.n};
    Tlwe output{param.n};
    symEncTlwe(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorusGeneral(mu, param.torusBase, LWE_Q) << endl;
    auto decPre = symDecTlweToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

    ScaledTlwe inputModN2{param.N * 2, param.n};
    Trlwe tv{param.k, param.N};
    Trlwe acc{param.k, param.N};
    std::vector<Trlwe8> accCRT(param.d, Trlwe8{param.k, param.N});
    Tlwe tmp{ksKey.nCurrKey};
//    COUNT_TIME("rescaleTlweFromTorus32", rescaleTlweFromTorus32(inputModN2, input);) // rescale to mod 2N
//    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(tv, v, inputModN2);) // tv = (X^-b) * (0,...,0,v)
//    COUNT_TIME("trlweMCRTDecomp", trlweMCRTDecomp(accCRT, tv, param);)
//    COUNT_TIME("blindRotateApproxCRTNtt", blindRotateApproxCRTNtt(accCRT, bsKeyCRT, inputModN2, param);)
//    COUNT_TIME("trlweMCRTToCRT", trlweMCRTToCRT(accCRT, param);)
//    COUNT_TIME("trlweCRTRecomp", trlweCRTRecomp(acc, accCRT, param);)
//    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, acc, param.driftPhase);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
//    COUNT_TIME("tlweKeySwitch", tlweKeySwitch(output, ksKey, tmp, param);)

    COUNT_TIME("functionalBootstrappingCrt", functionalBootstrappingCrt(output, input, bsKeyCRT, ksKey, v, param);)

    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;

    return 0;
}