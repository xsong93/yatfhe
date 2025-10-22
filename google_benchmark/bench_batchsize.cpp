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
        param.l = param.lApprox;
        initYatfhe(param);

        // key gen
        TlweKey tlweKey{param.n, param.lweStdDev};
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        bskMP = BootstrappingKeyMP{param};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
        TorusPolynomial v {param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);
        bskPre = BootstrappingKeyMPPreRot{param};
        genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorus32(pt, param.torusBase);
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
    BootstrappingKeyMP bskMP;
    BootstrappingKeyMPPreRot bskPre;
    Trlwe acc;
    Trlwe out;
    ScaledTlwe sTlwe;
};

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_MULTITHREAD)(benchmark::State& state) {
    const int batchSize = state.range(0);
    const int tasksPerThread = state.range(1);
    auto localParam = param;
    localParam.batchSize = batchSize;
    localParam.tasksPerThread = tasksPerThread;

    for (auto _ : state) {
        blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, localParam);
        benchmark::DoNotOptimize(out);
    }

    state.counters["batchSize"] = batchSize;
    state.counters["tasksPerThread"] = tasksPerThread;
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->ArgsProduct({
                          benchmark::CreateDenseRange(1, 20, 1),  // batchSize
                          benchmark::CreateDenseRange(1, 20, 1)   // tasksPerThread
    })
    ->ArgNames({"batchSize", "tasksPerThread"})
    ->MeasureProcessCPUTime()
    ->UseRealTime();

BENCHMARK_MAIN();