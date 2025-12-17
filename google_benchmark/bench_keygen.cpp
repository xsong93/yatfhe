#include <benchmark/benchmark.h>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

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
        state.PauseTiming();
        string file;
        file.append("BSK_GINX_")
            .append(to_string(state.iterations()))
            .append(".bin");
        serializeBskMP(bskMP, file);
        state.ResumeTiming();
    }
}

BENCHMARK_DEFINE_F(KeyGenBenchmark, LAZY_PIPE)(benchmark::State& state) {
    for (auto _ : state) {
        BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
        symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trgswKey.trlweKey, trgswKey.trlweKey.s, true, param);
        genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, v, param);
        state.PauseTiming();
        string file;
        file.append("BSK_LAZY_")
            .append(to_string(state.iterations()))
            .append(".bin");
        serializeBskLazyPipeAlt(bskMPLazyPipeAlt, file);
        state.ResumeTiming();
    }
}

BENCHMARK_REGISTER_F(KeyGenBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(KeyGenBenchmark, LAZY_PIPE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();

BENCHMARK_MAIN();