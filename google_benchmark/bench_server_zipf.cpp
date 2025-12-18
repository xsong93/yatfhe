#include <chrono>
#include <benchmark/benchmark.h>
#include <nlohmann/json.hpp>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/cache_manager.h"
#include "yautil/cache_workload_generator.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

YatfheParameters param{};
SimpleCacheManager cache(param.n, 5);
CacheWorkloadGenerator workload(0.8);
auto accessPattern = workload.generateAccessPattern(1000);

class ZipfBenchmark : public benchmark::Fixture {
public:
    ZipfBenchmark()
        : v{param.N},
        tlweKey{param.n, param.lweStdDev},
        trgswKey{param},
        dummyKsKey{param},
        input{param.n},
        acc{param},
        sTlwe{param.N * 2, param.n}
    {}

    void SetUp(const benchmark::State& state) override {
        initYatfhe(param);

        // key gen
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        generateTestPolynomial(v, param.torusBase, 2 * param.N);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(dummyKsKey, trlweKey, tlweKsKey, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        symEncTlwe(input, mu, tlweKey);
    }

protected:
    TorusPolynomial v;
    TlweKey tlweKey;
    TrgswKey trgswKey;
    TlweKeySwitchingKey dummyKsKey;
    Tlwe input;
    Trlwe acc;
    ScaledTlwe sTlwe;
};

BENCHMARK_DEFINE_F(ZipfBenchmark, GINX)(benchmark::State& state) {
    int request = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        clearFileCache();
        state.ResumeTiming();
        // rescaleTlweToNewMod(sTlwe, input);
        // genNoiselessTrlweSample(acc, v, sTlwe);
        auto& bskServer = cache.getGinxKey(accessPattern[request]);
        // blindRotateJP22Ntt(acc, bskServer.bskDft, sTlwe, param);
        // Tlwe tmp{dummyKsKey.nCurrKey}, output{param.n};
        // extractTlweFromTrlwe(tmp, acc, param.driftPhase);
        // switchKeyForTlwe(output, dummyKsKey, tmp, param);

        benchmark::DoNotOptimize(bskServer);
    }
}

BENCHMARK_DEFINE_F(ZipfBenchmark, LAZY)(benchmark::State& state) {
    int request = state.range(0);
    for (auto _ : state) {
        state.PauseTiming();
        clearFileCache();
        state.ResumeTiming();
        rescaleTlweToNewMod(sTlwe, input);
        auto* bskServer = cache.getLazyKeySimple(accessPattern[request]);
        Trlwe out{param};
        if (bskServer != nullptr) {
            blindRotateLazyPipeAltNtt(out, bskServer->bskFirst, bskServer->bskPrime,bskServer->s2Dft, sTlwe, v, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            std::string file = DiskReader::generateLazyKeyFilename(accessPattern[request]);
            blindRotateLazyPipeAltInitNtt(out, bsk.bskFirst, bsk.bskPrime,bsk.s2Dft, sTlwe,
                v, file, param);
            cache.putLazyKey(accessPattern[request], std::move(bsk));
        }
        Tlwe tmp{dummyKsKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, out, param.driftPhase);
        switchKeyForTlwe(output, dummyKsKey, tmp, param);

        benchmark::DoNotOptimize(output);
    }
}

BENCHMARK_REGISTER_F(ZipfBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->ArgsProduct({
        benchmark::CreateDenseRange(1, 100, 1)
    })
    ->ArgNames({"request"})
    ->UseRealTime();

BENCHMARK_REGISTER_F(ZipfBenchmark, LAZY)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->ArgsProduct({
        benchmark::CreateDenseRange(1, 100, 1)
    })
    ->ArgNames({"request"})
    ->UseRealTime();

BENCHMARK_MAIN();