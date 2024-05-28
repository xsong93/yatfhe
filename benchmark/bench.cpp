#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"
#include "yautil/control_helper.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yatfhe/keyswitching.h"

int main(int argc, char **argv) {
    YatfheParameters param {};
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);


    TlweKey tlweKey {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    COUNT_TIME("lweKeyGen", lweKeyGen(tlweKey);)
    COUNT_TIME("trlweKeyGen", trlweKeyGen(trlweKey);)
    COUNT_TIME("bootstrappingKeyGen", bootstrappingKeyGen(bsKey, param, trgswKey, tlweKey);)
    COUNT_TIME("tlweKeySwitchingKeyGen", tlweKeySwitchingKeyGen(ksKey, trlweKey, tlweKey, param);)

    Integer plain = 2;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorus32(mu, param.torusBase) << endl;
    auto decPre = symDecTlweSampleToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

//    COUNT_TIME("trgswFunctionalBootstrapping", trgswFunctionalBootstrappingNtt(output, input, bsKey, ksKey, v, param);)
    ScaledTlwe inputModN2 {param.N * 2, param.n};
    Trlwe accum {param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    COUNT_TIME("rescaleTlweFromTorus32", rescaleTlweFromTorus32(inputModN2, input);) // rescale to mod 2N
    COUNT_TIME("genNoiselessTrlweSample", genNoiselessTrlweSample(accum, v, inputModN2);) // accum = (X^-b) * (0,...,0,v)
    COUNT_TIME("blindRotateNtt", blindRotateNtt(accum, bsKey, inputModN2, param);)
    COUNT_TIME("extractTlweFromTrlwe", extractTlweFromTrlwe(tmp, accum, 0);) // tmp = (a', b0), a' = ((a1)0, -(a1)N-1, ... , -(a1)1, ..., ..., (ak)0, -(ak)N-1, ... , -(ak)1)
    COUNT_TIME("tlweKeySwitch", tlweKeySwitch(output, ksKey, tmp, param);)

    auto decAft = symDecTlweSampleToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;

    return 0;
}
