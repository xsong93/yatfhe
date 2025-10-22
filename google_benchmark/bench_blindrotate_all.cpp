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
        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);
        bskPre = BootstrappingKeyMPPreRot{param};
        genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);

        bskMPOpt = BootstrappingKeyMPOpt{param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);

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
    BootstrappingKeyMPOpt bskMPOpt;
    TorusPolynomial v;
    Trlwe acc;
    Trlwe out;
    ScaledTlwe sTlwe;
};

//BENCHMARK_DEFINE_F(BlindRotateBenchmark, MP21_SINGLETHREAD)(benchmark::State& state) {
//    for (auto _ : state) {
//        blindRotateMP21Ntt(acc, bskMP.bskDft, sTlwe, param);
//    }
//}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, JP22_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateJP22Ntt(acc, bskMP.bskDft, sTlwe, param);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateOptNtt(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, v, param);
    }
}

//BENCHMARK_DEFINE_F(BlindRotateBenchmark, JP22_MULTITHREAD)(benchmark::State& state) {
//    param.batchSize = state.range(0);
//    param.tasksPerThread = state.range(1);
//    for (auto _ : state) {
//        blindRotateJP22NttMT(acc, bskMP.bskDft, sTlwe, param);
//    }
//}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_MULTITHREAD)(benchmark::State& state) {
    param.batchSize = state.range(0);
    param.tasksPerThread = state.range(1);
    for (auto _ : state) {
        blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);
    }
}

//BENCHMARK_REGISTER_F(BlindRotateBenchmark, MP21_SINGLETHREAD)
//    ->Unit(benchmark::kMicrosecond)
//    ->Iterations(500)
//    ->MeasureProcessCPUTime()
//    ->UseRealTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, JP22_SINGLETHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
    ->MeasureProcessCPUTime()
    ->UseRealTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_SINGLETHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
    ->MeasureProcessCPUTime()
    ->UseRealTime();
//BENCHMARK_REGISTER_F(BlindRotateBenchmark, JP22_MULTITHREAD)
//    ->Unit(benchmark::kMicrosecond)
//    ->Iterations(500)
//        ->ArgsProduct({
//                              benchmark::CreateDenseRange(1, 20, 1),  // batchSize
//                              benchmark::CreateDenseRange(1, 20, 1)   // tasksPerThread
//                      })
//        ->ArgNames({"batchSize", "tasksPerThread"})
//    ->MeasureProcessCPUTime()
//    ->UseRealTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_MULTITHREAD)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
        ->ArgsProduct({
                              benchmark::CreateDenseRange(1, 20, 1),  // batchSize
                              benchmark::CreateDenseRange(1, 20, 1)   // tasksPerThread
                      })
        ->ArgNames({"batchSize", "tasksPerThread"})
    ->MeasureProcessCPUTime()
    ->UseRealTime();

BENCHMARK_MAIN();