#include <sys/time.h>
#include <iostream>
#include "yatfhe.h"
#include "time_counter.h"
#include "parameters.h"
#include "benchmark.h"

// Benchmark params
#define _EXECS 100

// Benchmark functions (some depend on others)
#define BENCH_TRGSW_BOOTSTRAP

int main(int argc, char **argv) {
    auto* timer = TimeCounter::init();
    COUNT_TIME("test", timer, std::cout << n << std::endl;)
    COUNT_TIME("test", timer, std::cout << N << std::endl;)
    COUNT_TIME("test", timer, std::cout << k << std::endl;)
    free(timer);
//    TLWE_Key key_tlwe = tlwe_new_binary_key(n, lwe_std_dev);
//    TLWE_Key key_tlwe_out = tlwe_new_binary_key(N, rlwe_std_dev);
//    TRLWE_Key key_trlwe = trlwe_new_binary_key(N, k, rlwe_std_dev);
//    trlwe_extract_tlwe_key(key_tlwe_out, key_trlwe);
//    TRGSW_Key trgsw_key = trgsw_new_key(key_trlwe, l, Bg_bit);
//    TLWE_KS_Key tlwe_ksk = tlwe_new_KS_key(key_tlwe, key_tlwe_out, t, base_bit);
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
//    Bootstrap_Key bk_key = new_bootstrap_key(trgsw_key, key_tlwe, 1);
//    Bootstrap_GA_Key bk_ga_key = new_bootstrap_key_ga(trgsw_key, key_tlwe);
//
//#ifdef BENCH_TRGSW_BOOTSTRAP
//    TRGSW_DFT c_trgsw = trgsw_alloc_new_DFT_sample(l, Bg_bit, k, N);
//    BENCHMARK("TRGSW_BS_P1", _EXECS, "TRGSW Functional Bootstrap PHASE 1",
//    functional_bootstrap_trgsw_phase1(c_trgsw, sel, bk_key, 4);
//    );
//
//    BENCHMARK("TRGSW_BS_P2", _EXECS, "TRGSW Functional Bootstrap PHASE 2",
//    functional_bootstrap_trgsw_phase2(c[4], c_trgsw, lut_c);
//    );
//#endif

    return 0;
}
