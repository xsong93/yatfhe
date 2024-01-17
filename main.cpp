#include <iostream>
#include "tlwe.h"
#include "trlwe.h"
#include "trgsw.h"
#include "bootstrap.h"
#include "time_counter.h"
#include "yatfhe_parameters.h"
#include "control_helper.h"

int main(int argc, char **argv) {
    auto* timer = TimeCounter::init();
    const YatfheParameters yatfheParameters {};
    COUNT_TIME("test n", timer, std::cout << yatfheParameters.n << std::endl;)
    COUNT_TIME("test N", timer, std::cout << yatfheParameters.N << std::endl;)
    COUNT_TIME("test k", timer, std::cout << yatfheParameters.k << std::endl;)

    TlweKey tlweKey {};
    newBinaryTlweKey(tlweKey, yatfheParameters);
    TlweKey keyTlweOut {};
    newBinaryTlweKey(keyTlweOut, yatfheParameters);
    TrlweKey trlweKey {};
    newBinaryTrlweKey(trlweKey, yatfheParameters);
//    trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
    TrgswKey trgswKey {trlweKey, yatfheParameters.l, yatfheParameters.bgBit};

//    TLWE_KS_Key tlwe_ksk = tlwe_new_KS_key(key_tlwe, key_tlwe_out, t, baseBit);
//
//    auto * input = static_cast<Torus *>(safe_aligned_malloc(sizeof(Torus) * (_EXECS * 4 + 1)));
//    generate_random_bytes(sizeof(Torus)*(_EXECS*4 + 1), (uint8_t *) input);
//    TLWE c[_EXECS*4 + 1];
//    for (size_t i = 0; i < _EXECS*4 + 1; i++) c[i] = tlwe_new_sample(input[i], key_tlwe_out);
//
//    TLWE sel = tlwe_new_sample(double2torus(1./8), key_tlwe);
//    TorusPolynomial poly_res = polynomial_new_torus_polynomial(N);
//    for (size_t i = 0; i < N; i++) poly_res->coeffs[i] = input[i / (N / 4)];
//
//    TRLWE lut_c = trlwe_new_noiseless_trivial_sample(poly_res, k, N);

    BootstrappingKey bsKey {};
    newBootstrappingKey(bsKey, yatfheParameters, trgswKey, tlweKey);
//    Bootstrap_GA_Key bk_ga_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);
//
//#ifdef BENCH_TRGSW_BOOTSTRAP
    TrgswDft newDftSample {};
    initTrgswDftSample(newDftSample, yatfheParameters);
////    BENCHMARK("TRGSW_BS_P1", _EXECS, "TRGSW Functional Bootstrap PHASE 1",
//    functional_bootstrap_trgsw_phase1(newDftSample, sel, bsKey, 4);
////    );
//
////    BENCHMARK("TRGSW_BS_P2", _EXECS, "TRGSW Functional Bootstrap PHASE 2",
//    functional_bootstrap_trgsw_phase2(c[4], newDftSample, lut_c);
////    );
////#endif

    std::cout <<1;
//    deleteTrgswKey(trgswKey);
    delete timer;
    return 0;
}
