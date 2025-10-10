#include <benchmark/benchmark.h>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric_functions.h"
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
        param.batchSize = 40;
        param.l = 2;
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
        rescaleTlweFromTorus32(sTlwe, input);
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

BENCHMARK_DEFINE_F(BlindRotateBenchmark, MP21_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateMP21Ntt(acc, bskMP.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(acc);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, JP22_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateJP22Ntt(acc, bskMP.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(acc);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateWithPreRotNtt(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(out);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, JP22_MULTITHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateJP22NttMT(acc, bskMP.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(acc);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_MULTITHREAD)(benchmark::State& state) {
    const int batchSize = state.range(0);
    auto localParam = param;
    localParam.batchSize = batchSize;

    for (auto _ : state) {
        blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, localParam);
        benchmark::DoNotOptimize(out);
    }

    state.counters["batchSize"] = batchSize;
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, MP21_SINGLETHREAD)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, JP22_SINGLETHREAD)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_SINGLETHREAD)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, JP22_MULTITHREAD)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->DenseRange(1, 40, 1)
    ->ArgNames({"batchSize"});

BENCHMARK_MAIN();