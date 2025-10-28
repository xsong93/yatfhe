#include <thread>
#include <benchmark/benchmark.h>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"

class BlindRotateBenchmark : public benchmark::Fixture {
public:
    BlindRotateBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        param = YatfheParameters{};
        param.N = 512;
        initYatfhe(param);

        // key gen
        TlweKey tlweKey{param.n, param.lweStdDev};
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        TlweKeySwitchingKey ksKey{param};
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
        s2Dft = TrlevDft{param, param.l};
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        bskMPLazy = BootstrappingKeyMPOpt{param, param.lApprox, true};
        genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        sTlwe = ScaledTlwe{param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        acc = Trlwe{param.k, param.N};
        genNoiselessTrlweSample(acc, v, sTlwe);

        out = Trlwe{param.k, param.N};
    }

protected:
    YatfheParameters param;
    BootstrappingKeyMPOpt bskMPLazy;
    Trlwe acc;
    Trlwe out;
    TrlevDft s2Dft;
    ScaledTlwe sTlwe;
    TorusPolynomial v;
};

BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_MULTITHREAD)(benchmark::State& state) {
    const int batchSize = state.range(0);
    const int tasksPerThread = state.range(1);
    auto localParam = param;
    localParam.batchSize = batchSize;
    localParam.tasksPerThread = tasksPerThread;

    for (auto _ : state) {
        blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft, bskMPLazy.initialized, sTlwe, v, s2Dft, localParam);
        benchmark::DoNotOptimize(out);
    }

    state.counters["batchSize"] = batchSize;
    state.counters["tasksPerThread"] = tasksPerThread;
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->ArgsProduct({
                          benchmark::CreateDenseRange(1, 20, 1),  // batchSize
                          benchmark::CreateDenseRange(1, 20, 1)   // tasksPerThread
    })
    ->ArgNames({"batchSize", "tasksPerThread"})
    ->MeasureProcessCPUTime()
    ->UseRealTime();

BENCHMARK_MAIN();