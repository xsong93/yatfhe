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
    TlweKeySwitchingKey ksKey {param.N * param.k, param.n, param.ksLevel};
    lweKeyGen(tlweKey);
//    lweKeyGen(keyTlweOut, param.n);
    trlweKeyGen(trlweKey);
    bootstrappingKeyGen(bsKey, param, trgswKey, tlweKey);
    tlweKeySwitchingKeyGen(ksKey, trlweKey, tlweKey, param);

    Torus mu = doubleToTorus32(1.0 / param.torusBase);
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);
    Tlwe input {param.n};
    Tlwe output {param.n};
    symEncTlweSample(input, mu, tlweKey);

    cout <<"msg:"<<torus32ToDouble(mu)<<endl;
    auto decPre = symDecTlweSample(input, tlweKey, param.torusBase);
    cout <<"decPre:"<< decPre <<endl;

    trgswFunctionalBootstrapping(output, input, bsKey, ksKey, v, param);
//    printTlweAB(input, "input boot");
//    printTlweAB(output, "output boot");
    double decAft = symDecTlweSample(output, tlweKey, param.torusBase);
    cout <<"decAft:"<<decAft<<endl;


//    trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
//    TLWE_KS_Key tlwe_ksk = tlwe_new_KS_key(key_tlwe, key_tlwe_out, t, baseBit);
//    auto * input = static_cast<Torus *>(safe_aligned_malloc(sizeof(Torus) * (_EXECS * 4 + 1)));
//    generate_random_bytes(sizeof(Torus)*(_EXECS*4 + 1), (uint8_t *) input);
//    TLWE c[_EXECS*4 + 1];
//    for (size_t i = 0; i < _EXECS*4 + 1; i++) c[i] = tlwe_new_sample(input[i], key_tlwe_out);
//
//    TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
//    for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = input[i / (N / 4)];
//
//    TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);
//    Bootstrap_GA_Key bk_ga_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);
//
//#ifdef BENCH_TRGSW_BOOTSTRAP

//    initTrgswDftSample(output, param);
////    BENCHMARK("TRGSW_BS_P1", _EXECS, "TRGSW Functional Bootstrap PHASE 1",

////    );
//
////    BENCHMARK("TRGSW_BS_P2", _EXECS, "TRGSW Functional Bootstrap PHASE 2",
//    functional_bootstrap_trgsw_phase2(c[4], output, lut_c);
////    );
////#endif

    return 0;
}
