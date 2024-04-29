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

int main(int argc, char **argv) {
//    TimeCounter timer {};
    const YatfheParameters param {};
    COUNT_TIME("init timer", std::cout << std::endl;)

    TlweKey tlweKey {param.n, param.lweStdDev};
//    TlweKey keyTlweOut {param.n, param.lweStdDev};
    TrgswKey trgswKey {param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey {param};
    TlweKeySwitchingKey ksKey {param};
    lweKeyGen(tlweKey);
//    lweKeyGen(keyTlweOut, param.n);
    trlweKeyGen(trlweKey);
    bootstrappingKeyGen(bsKey, param, trgswKey, tlweKey);
    tlweKeySwitchingKeyGen(ksKey, trlweKey, tlweKey, param);

    Integer plain = 2;
    Torus mu = modSwitchToTorus32(plain, param.torusBase);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    std::vector<Integer> t(v.N);
    for (int i = 0; i < v.N; i++) {
        t[i] = modSwitchFromTorus32(v.coeffs[i], param.torusBase);
    }

    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout << "msg: " << modSwitchFromTorus32(mu, param.torusBase) << endl;
    auto decPre = symDecTlweSampleToInt(input, tlweKey, param.torusBase);
    cout << "decPre: " << decPre << endl;

    trgswFunctionalBootstrapping(output, input, bsKey, ksKey, v, param);

    auto decAft = symDecTlweSampleToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;

    return 0;
}
