#include <benchmark/benchmark.h>
#include <filesystem>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

class BlindRotateBenchmark : public benchmark::Fixture {
public:
    BlindRotateBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        initYatfhe(param);

        // key gen
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        rescaleTlweToNewMod(sTlwe, input);
        genNoiselessTrlweSample(acc, v, sTlwe);
        encryptTrgswMPNtt(one, 1, trgswKey, 0, param);
    }

protected:
    YatfheParameters param{};
    TorusPolynomial v{param.N};
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    Trlwe acc{param};
    Trlwe out{param};
    ScaledTlwe sTlwe{param.N * 2, param.n};
    TrlevDft s2Dft{param, param.l};
    TrgswMPDft one{param};
};

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX)(benchmark::State& state) {
    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
    for (auto _ : state) {
        blindRotateJP22Ntt(acc, bskMP.bskDft, sTlwe, param);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX_OPT)(benchmark::State& state) {
    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    for (auto _ : state) {
        blindRotateOptNtt(out, bskMPOpt.bskFirst, bskMPOpt.bskDft, sTlwe, v, param);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_SINGLETHREAD)(benchmark::State& state) {
    BootstrappingKeyMPLazy bskMPLazyOpt{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);
    auto localP = param;
    localP.batchSize = 1;
    localP.tasksPerThread = 1;
    for (auto _ : state) {
        blindRotateLazyMTNtt(out, bskMPLazyOpt.bskFirst, bskMPLazyOpt.bskTrim, bskMPLazyOpt.bskDecompA,
                             sTlwe, v, s2Dft, localP);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_NAIVE_MULTITHREAD)(benchmark::State& state) {
    BootstrappingKeyMPLazy bskMPLazyOpt{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);
    auto localP = param;
    localP.batchSize = 3;
    localP.tasksPerThread = 1;
    for (auto _ : state) {
        blindRotateLazyMTNtt(out, bskMPLazyOpt.bskFirst, bskMPLazyOpt.bskTrim, bskMPLazyOpt.bskDecompA,
                             sTlwe, v, s2Dft, localP);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_PIPELINE)(benchmark::State& state) {
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    for (auto _ : state) {
        blindRotateLazyPipeNtt(out, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft,
                               bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one, param);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, LAZY_PIPELINE_S)(benchmark::State& state) {
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    for (auto _ : state) {
        blindRotateLazyPipeSerializationNtt(out, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft,
                               bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one, "LAZY_PIPELINE_S.bin",
                               false, param);
    }
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_SINGLETHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_NAIVE_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_PIPELINE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, LAZY_PIPELINE_S)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime()
    ->MeasureProcessCPUTime();

BENCHMARK_MAIN();