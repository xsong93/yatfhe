#include <benchmark/benchmark.h>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"

class KeyGenBenchmark : public benchmark::Fixture {
public:
    KeyGenBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        initYatfhe(param);

        // key gen
        TrlweKey& trlweKey = trgswKey.trlweKey;
        TlweKeySwitchingKey ksKey{param};
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
        generateTestPolynomial(v, param.torusBase, 2 * param.N);
    }

protected:
    YatfheParameters param{};
    TrgswKey trgswKey{param};
    TlweKey tlweKey{param.n, param.lweStdDev};
    TorusPolynomial v{param.N};
};

BENCHMARK_DEFINE_F(KeyGenBenchmark, GINX)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMP bskMP{param, param.lApprox};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
    }
}

BENCHMARK_DEFINE_F(KeyGenBenchmark, GINX_OPT)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMPOpt bskMPOpt {param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    }
}

BENCHMARK_DEFINE_F(KeyGenBenchmark, OURS)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
        genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    }
}

BENCHMARK_REGISTER_F(KeyGenBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(KeyGenBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(KeyGenBenchmark, OURS)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();

BENCHMARK_MAIN();