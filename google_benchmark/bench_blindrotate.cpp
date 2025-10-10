#include <benchmark/benchmark.h>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric_functions.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

// 继承 benchmark::Fixture
class BlindRotateBenchmark : public benchmark::Fixture {
public:
    BlindRotateBenchmark() = default;

    // 在每个测试用例运行前执行，用于初始化资源
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

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateMPNtt(acc, bskMP.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(acc);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_SINGLETHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateWithPreRotNtt(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(out);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX_MULTITHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateMPNttMT(acc, bskMP.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(out);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_MULTITHREAD)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);
        benchmark::DoNotOptimize(out);
    }
}

// 注册测试
BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_SINGLETHREAD);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_SINGLETHREAD);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_MULTITHREAD);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_MULTITHREAD);

// 使用库提供的main函数
BENCHMARK_MAIN();