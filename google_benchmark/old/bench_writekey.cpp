#include <benchmark/benchmark.h>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

void serializeBskMPApp(const BootstrappingKeyMP& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::app);
    for (auto i = 0; i < t.n; i++) {
        serialize(t.bskDft[i][0], os);
    }
    writePOD(os, t.n);
    writePOD(os, t.group);
    writePOD(os, t.isHalf);
    os.close();
}

void serializeBskMPOptApp(const BootstrappingKeyMPOpt& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::app);
    serialize(t.bskFirst[0], os);
    for (auto i = 0; i < t.n - 1; i++) {
        serialize(t.bskDft[i][0], os);
    }
    writePOD(os, t.n);
    writePOD(os, t.group);
    os.close();
}

void serializeBskLazyPipeApp(BootstrappingKeyMPLazyPipe& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::app);

    // Step 1: Write bskFirst[0] and bskDft[0][0] (written outside loop)
    serialize(t.bskFirst[0], os);
    serialize(t.bskDft[0][0], os);

    // Step 2: Write alternating pattern: bskDft[i+1][0], bskDecompA[i][0]
    for (int i = 0; i < t.n - 2; ++i) {
        // Write bskDft[i+1][0]
        serialize(t.bskDft[i + 1][0], os);

        // Write bskDecompA[i][0] (if applicable)
        if (i < t.n - 3) {
            serializeNestedVector(t.bskDecompA[i][0], os);
        }
    }
    os.close();
}

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
    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
    for (auto _ : state) {
        serializeBskMPApp(bskMP, "BSK_GINX_W_BENCH.bin");
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, GINX_OPT)(benchmark::State& state) {
    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    for (auto _ : state) {
        serializeBskMPOptApp(bskMPOpt, "BSK_GINX_OPT_W_BENCH.bin");
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, LAZY_PIPE)(benchmark::State& state) {
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    for (auto _ : state) {
        serializeBskLazyPipeApp(bskMPLazyPipe, "BSK_PIPE_W_BENCH.bin");
    }
}

BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, LAZY_PIPE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(10)
    ->UseRealTime();

BENCHMARK_MAIN();