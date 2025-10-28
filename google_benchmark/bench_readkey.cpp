#include <benchmark/benchmark.h>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

class ReadKeyBenchmark : public benchmark::Fixture {
public:
    ReadKeyBenchmark() = default;

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

BENCHMARK_DEFINE_F(ReadKeyBenchmark, GINX)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMP bskMPServer;
        deserializeBskMP(bskMPServer, "bsk_serialized_GINX.bin", param.n);
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, GINX_OPT)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMPOpt bskMPOpt {param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, OURS)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
        genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    }
}

BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, OURS)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime()
    ->MeasureProcessCPUTime();

BENCHMARK_MAIN();